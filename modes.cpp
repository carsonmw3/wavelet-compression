#include "modes.h"

#include "preprocess.h"
#include "box-structs.h"
#include "readandwrite.h"
#include "compressor.h"
#include "decompressor.h"
#include "calc-loss.h"
#include "writeplotfile.h"
#include "iterator.h"
#include "tmpdir.h"

#include <numeric>
#include <spdlog/spdlog.h>
#include <chrono>
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>


int compress(const Config& cfg) {

    std::vector<std::string> files  = format_files(cfg.data_dir, cfg.min_time, cfg.max_time);
    std::vector<int>         levels = format_levels(cfg.min_level, cfg.max_level);

    // total number of timesteps, levels, and components, for use in iteration
    int num_times = files.size();
    int num_levels = levels.size();
    int num_components = cfg.components.size();

    spdlog::info("Processing data...");
    auto start = std::chrono::high_resolution_clock::now();

    RunInfo runinfo;
    runinfo.files      = files;
    runinfo.min_level  = cfg.min_level;
    runinfo.max_level  = cfg.max_level;
    runinfo.components = cfg.components;

    // process data for compression
    AllData data = preprocess_data(files,
                                   cfg.components,
                                   levels);

    auto& boxes       = data.boxes;
    auto& locations   = data.locations;
    auto& dimensions  = data.dimensions;
    auto& box_counts  = data.box_counts;
    auto& min_values  = data.min_values;
    auto& max_values  = data.max_values;
    auto& amrexinfo   = data.amrexinfo;

    runinfo.comp_idxs = data.comp_idxs;

    AMRIterator iterator(num_times, num_levels, box_counts, num_components);

    // Verify access to compresseddir
    if (!cfg.compressed_dir.empty() && !std::filesystem::exists(cfg.compressed_dir)) {
        std::error_code ec;
        std::filesystem::create_directories(cfg.compressed_dir, ec);
        if (ec) {
            spdlog::error("Failed to create compressed directory {}: {}",
                          cfg.compressed_dir, ec.message());
        }
    }

    // write runinfo, location, dimension, box count, amrex data
    write_runinfo(runinfo,
                  cfg.compressed_dir,
                  "runinfo.raw");
    write_loc_dim_to_bin(locations,
                         cfg.compressed_dir,
                         "locations.raw",
                         iterator);
    write_loc_dim_to_bin(dimensions,
                         cfg.compressed_dir,
                         "dimensions.raw",
                         iterator);
    write_box_counts(box_counts,
                     cfg.compressed_dir,
                     "boxcounts.raw",
                     num_times,
                     num_levels);
    write_amrexinfo(amrexinfo,
                    cfg.compressed_dir,
                    "amrexinfo.raw");

    auto end      = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end - start).count();
    spdlog::info(
        "Successfully processed data in {} seconds. Beginning compression...",
        duration);

    auto start1 = std::chrono::high_resolution_clock::now();

    // compress
    iterator.iterate([&](int t, int lev, int box_idx) {
        multiBox3D& current_box = boxes[t][lev][box_idx];
        compress(current_box, runinfo.comp_idxs, cfg.keep, t, lev, box_idx, cfg.compressed_dir);
    });

    auto end1      = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration<double>(end1 - start1).count();
    spdlog::info(
        "Compression completed in {} seconds.",
        duration1);

    return 0;
}


int decompress(const Config& cfg) {

    RunInfo runinfo = read_runinfo(cfg.compressed_dir, "runinfo.raw");

    std::vector<int> levels = format_levels(runinfo.min_level, runinfo.max_level);

    spdlog::info("Decompressing data between timestep {} and {}, level {} and {}, for {} components",
                 runinfo.files[0], runinfo.files[runinfo.files.size() - 1], runinfo.min_level, runinfo.max_level,
                 runinfo.components.size());

    // total number of timesteps, levels, and components, for use in iteration
    int num_times = runinfo.files.size();
    int num_levels = levels.size();
    int num_components = runinfo.components.size();

    spdlog::info("Beginning decompression...");
    auto start2 = std::chrono::high_resolution_clock::now();

    // get number of boxes at each time and level
    auto counts_read = read_box_counts(cfg.compressed_dir,
                                       "boxcounts.raw",
                                       num_times,
                                       num_levels);

    AMRIterator iterator(num_times, num_levels, counts_read, num_components);

    // for storage of regenerated boxes
    std::vector<std::vector<std::vector<multiBox3D>>> regen_boxes;

    // reshape according to number of timesteps and levels
    regen_boxes.resize(num_times);
    for (int t = 0; t < num_times; ++t) {
        regen_boxes[t].resize(num_levels);
    }

    // decompress
    iterator.iterate([&](int t, int lev, int box_idx) {
        if (regen_boxes[t][lev].size() == 0)
            regen_boxes[t][lev].resize(counts_read[t][lev]);

        multiBox3D multibox;
        for (int component : runinfo.comp_idxs) {
            std::string file_path = cfg.compressed_dir + "compressed-wavelet-" +
                                    std::to_string(t) + "-" +
                                    std::to_string(lev) + "-" +
                                    std::to_string(component) + "-" +
                                    std::to_string(box_idx) + ".xz";
            Box3D curr_box = decompress(file_path, t, lev, component, box_idx);
            multibox.push_back(std::move(curr_box));
        }
        regen_boxes[t][lev][box_idx] = std::move(multibox);
    });

    auto end2      = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration<double>(end2 - start2).count();
    spdlog::info("Decompression completed in {} seconds.",
                 duration2);

    // read info required for writing plotfiles, then write them
    AMReXInfo amrexinfo_read = read_amrex_info(cfg.compressed_dir,
                                               "amrexinfo.raw");

    LocDimData locs_read = read_loc_dim_from_bin(cfg.compressed_dir,
                                                 "locations.raw",
                                                 counts_read,
                                                 iterator,
                                                 num_times,
                                                 num_levels);

    LocDimData dims_read = read_loc_dim_from_bin(cfg.compressed_dir,
                                                 "dimensions.raw",
                                                 counts_read,
                                                 iterator,
                                                 num_times,
                                                 num_levels);

    write_plotfiles(std::move(regen_boxes),
                    locs_read,
                    dims_read,
                    runinfo.files,
                    num_levels,
                    num_components,
                    runinfo.components,
                    amrexinfo_read,
                    "../../regenerated-plotfiles/");

    spdlog::info("Sucessfully wrote plotfiles.");

    return 0;
}


// runs a quick compression/decompression with a limited amount of data to
// get an estimate for metrics like rmse, loss, and compressed size
int estimate(Config& cfg) {

    // only runs on one time/level, but for all components
    int num_times = 1;
    int num_levels = 1;
    int num_components = cfg.components.size();

    // storage for compressed files used in estimation
    TempDir scratch_dir;

    std::vector<std::string> files = format_files(cfg.data_dir, cfg.min_time, cfg.min_time);

    std::vector<int> levels;
    levels.push_back(cfg.min_level);

    AllData data = preprocess_data(files,
                                   cfg.components,
                                   levels);

    auto& boxes       = data.boxes;
    auto& box_counts  = data.box_counts;
    auto& min_values  = data.min_values;
    auto& max_values  = data.max_values;
    auto& amrexinfo   = data.amrexinfo;

    AMRIterator iterator(num_times, num_levels, box_counts, num_components);

    iterator.iterate([&](int t, int lev, int box_idx) {
        multiBox3D& current_box = boxes[t][lev][box_idx];
        compress(current_box, data.comp_idxs, cfg.keep, t, lev, box_idx, scratch_dir.path());
    });

    spdlog::info("Compression complete.");

    std::vector<std::vector<std::vector<multiBox3D>>> regen_boxes;

    regen_boxes.resize(num_times);
    for (int t = 0; t < num_times; ++t) {
        regen_boxes[t].resize(num_levels);
    }

    iterator.iterate([&](int t, int lev, int box_idx) {
        if (regen_boxes[t][lev].size() == 0)
            regen_boxes[t][lev].resize(box_counts[t][lev]);

        multiBox3D multibox;
        for (int component : data.comp_idxs) {
            std::filesystem::path file_path = scratch_dir.path() /
                                              ("compressed-wavelet-" + std::to_string(t) + "-" +
                                               std::to_string(lev) + "-" +
                                               std::to_string(component) + "-" +
                                               std::to_string(box_idx) + ".xz");
            Box3D curr_box = decompress(file_path, t, lev, component, box_idx);
            multibox.push_back(std::move(curr_box));
        }
        regen_boxes[t][lev][box_idx] = std::move(multibox);
    });

    spdlog::info("Decompression complete.");

    std::vector<std::vector<double>> all_rmses(num_components);

    iterator.iterate([&](int t, int l, int b) {
        const auto& actual = boxes[t][l][b];
        const auto& regen  = regen_boxes[t][l][b];

        std::vector<double> rmse = calc_rmse_per_box(actual, regen, num_components);

        for (int c = 0; c < num_components; ++c) {
            all_rmses[c].push_back(rmse[c]);
        }
    });

           // Compute mean RMSE and adjusted loss
    for (int c = 0; c < num_components; ++c) {
        double mean_rmse = std::accumulate(all_rmses[c].begin(), all_rmses[c].end(), 0.0) /
                           all_rmses[c].size();

        spdlog::info("Predicted RMSE, {} = {}", cfg.components[c], mean_rmse);

        double loss = calc_adj_loss(mean_rmse, max_values[c] - min_values[c]);
        spdlog::info("Predicted Adjusted loss, {} = {}", cfg.components[c], loss);
    }

           // Estimate compressed size
    std::ifstream x;

           // Get total number of components in original data
    std::string header = files[0] + "/Header";
    x.open(header.c_str(), std::ios::in);
    if (!x.is_open()) {
        spdlog::error("Failed to open header file: {}", files[0]);
    }

           // read in first line of header
    std::string str;
    x >> str;

           // read in number of components from header
    int nComp;
    x >> nComp;
    x.close();

           // calculate original data size
    std::string raw_path = files[0] + "/Level_" + std::to_string(levels[0]) + "/";
    double raw_size = calc_size(raw_path);

           // adjust based on number of components being compressed
    raw_size /= nComp;
    raw_size *= num_components;

           // calculate compressed size
    double compressed_size = calc_size(scratch_dir.path());

           // print as percentage of original size
    spdlog::info("Predicted compressed size: {}%", (compressed_size / raw_size) * 100);

    return 0;

}

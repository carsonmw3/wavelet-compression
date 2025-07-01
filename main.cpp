#include "argparse.h"
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

    // total number of timesteps, levels, and components, for use in iteration
    int num_times = cfg.max_time - cfg.min_time + 1;
    int num_levels = cfg.max_level - cfg.min_level + 1;
    int num_components = cfg.components.size();

    // vector of name of directory for each timestep
    std::vector<std::string> files;
    for (int t = cfg.min_time; t <= cfg.max_time; ++t)
        files.push_back(cfg.data_dir + "plt0" + std::to_string(t) + "00");

    // vector of levels
    std::vector<int> levels;
    for (int l = cfg.min_level; l <= cfg.max_level; ++l)
        levels.push_back(l);

    spdlog::info("Processing data...");
    auto start = std::chrono::high_resolution_clock::now();

    RunInfo runinfo;
    runinfo.min_time = cfg.min_time;
    runinfo.max_time = cfg.max_time;
    runinfo.min_level = cfg.min_level;
    runinfo.max_level = cfg.max_level;
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
    data.runinfo      = runinfo;

    AMRIterator iterator(num_times, num_levels, box_counts, num_components);

    // write runinfo, location, dimension, box count, amrex data
    write_runinfo(data.runinfo,
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
        compress(current_box, cfg.components, cfg.keep, t, lev, box_idx, cfg.compressed_dir);
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

    spdlog::info("Decompressing data between timestep {} and {}, level {} and {}, for {} components",
                 runinfo.min_time, runinfo.max_time, runinfo.min_level, runinfo.max_level,
                 runinfo.components.size());

    // total number of timesteps, levels, and components, for use in iteration
    int num_times = runinfo.max_time - runinfo.min_time + 1;
    int num_levels = runinfo.max_level - runinfo.min_level + 1;
    int num_components = runinfo.components.size();

    spdlog::info("Beginning decompression...");
    auto start2      = std::chrono::high_resolution_clock::now();

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
        for (int component : runinfo.components) {
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
    spdlog::info("Decompression completed in {} seconds. Writing plotfiles...",
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

    write_plotfiles(regen_boxes,
                    locs_read,
                    dims_read,
                    num_times,
                    num_levels,
                    num_components,
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

    std::vector<std::string> files;
    files.push_back(cfg.data_dir + "plt0" + std::to_string(cfg.min_time) + "00");

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
        compress(current_box, cfg.components, cfg.keep, t, lev, box_idx, scratch_dir.path());
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
        for (int component : cfg.components) {
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

        spdlog::info("Predicted RMSE, {} = {}", amrexinfo.comp_names[c], mean_rmse);

        double loss = calc_adj_loss(mean_rmse, max_values[c] - min_values[c]);
        spdlog::info("Predicted Adjusted loss, {} = {}", amrexinfo.comp_names[c], loss);
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


bool has_flag(int argc, char* argv[], const std::string& flag) {
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == flag) {
            return true;
        }
    }
    return false;
}


int main(int argc, char* argv[]) {

    amrex::Initialize(argc, argv);
    spdlog::set_level(spdlog::level::debug);

    if (has_flag(argc, argv, "-c")) {
        Config cfg = parse_config_compress();
        compress(cfg);
    } else if (has_flag(argc, argv, "-estimate")) {
        Config cfg = parse_config_compress();
        estimate(cfg);
    } else {
        Config cfg = parse_config_decompress();
        decompress(cfg);
    }

    amrex::Finalize();
    return 0;

}

    // std::vector<std::string> files;
    // for (int t = cfg.min_time; t < cfg.max_time+1; t++) {
    //     std::string filename = cfg.data_dir +"plt0" + std::to_string(t) + "00";
    //     files.push_back(filename);
    // }

    // std::vector<int> levels;
    // for (int l = cfg.min_level; l < cfg.max_level+1; l++) {
    //     levels.push_back(l);
    // }

    // spdlog::info("Processing data...");
    // auto start = std::chrono::high_resolution_clock::now();

    // int num_times      = cfg.max_time - cfg.min_time + 1;
    // int num_levels     = cfg.max_level - cfg.min_level + 1;
    // int num_components = cfg.components.size();

    // AllData data = preprocess_data(files,
    //                                cfg.components,
    //                                levels);

    // auto& boxes       = data.boxes;
    // auto& locations   = data.locations;
    // auto& dimensions  = data.dimensions;
    // auto& box_counts  = data.box_counts;
    // auto& min_values  = data.min_values;
    // auto& max_values  = data.max_values;
    // auto& amrexinfo   = data.amrexinfo;

    // AMRIterator iterator(num_times, num_levels, box_counts, num_components);

    // write_loc_dim_to_bin(locations,
    //                      cfg.compressed_dir,
    //                      "locations.raw",
    //                      iterator);
    // write_loc_dim_to_bin(dimensions,
    //                      cfg.compressed_dir,
    //                      "dimensions.raw",
    //                      iterator);
    // write_box_counts(box_counts,
    //                  cfg.compressed_dir,
    //                  "boxcounts.raw",
    //                  num_times,
    //                  num_levels);

    // write_amrexinfo(amrexinfo,
    //                 cfg.compressed_dir,
    //                 "amrexinfo.raw");

    // auto end      = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration<double>(end - start).count();
    // spdlog::info(
    //     "Successfully processed data in {} seconds. Beginning compression...",
    //     duration);

    // auto start1 = std::chrono::high_resolution_clock::now();

    // iterator.iterate([&](int t, int lev, int box_idx) {
    //     multiBox3D& current_box = boxes[t][lev][box_idx];
    //     compress(current_box, cfg.components, cfg.keep, t, lev, box_idx, cfg.compressed_dir);
    // });

    // auto end1      = std::chrono::high_resolution_clock::now();
    // auto duration1 = std::chrono::duration<double>(end1 - start1).count();
    // spdlog::info(
    //     "Compression completed in {} seconds. Beginning decompression...",
    //     duration1);

    // auto start2      = std::chrono::high_resolution_clock::now();
    // auto counts_read = read_box_counts(cfg.compressed_dir,
    //                                    "boxcounts.raw",
    //                                    num_times,
    //                                    num_levels);

    // std::vector<std::vector<std::vector<multiBox3D>>> regen_boxes;

    // regen_boxes.resize(num_times);
    // for (int t = 0; t < num_times; ++t) {
    //     regen_boxes[t].resize(num_levels);
    // }

    // iterator.iterate([&](int t, int lev, int box_idx) {
    //     if (regen_boxes[t][lev].size() == 0)
    //         regen_boxes[t][lev].resize(box_counts[t][lev]);

    //     multiBox3D multibox;
    //     for (int component : cfg.components) {
    //         std::string file_path = cfg.compressed_dir + "compressed-wavelet-" +
    //                                 std::to_string(t) + "-" +
    //                                 std::to_string(lev) + "-" +
    //                                 std::to_string(component) + "-" +
    //                                 std::to_string(box_idx) + ".xz";
    //         Box3D curr_box = decompress(file_path, t, lev, component, box_idx);
    //         multibox.push_back(std::move(curr_box));
    //     }
    //     regen_boxes[t][lev][box_idx] = std::move(multibox);
    // });

    // auto end2      = std::chrono::high_resolution_clock::now();
    // auto duration2 = std::chrono::duration<double>(end2 - start2).count();
    // spdlog::info("Decompression completed in {} seconds. Calculating loss...",
    //              duration2);

    // AMReXInfo amrexinfo_read = read_amrex_info(cfg.compressed_dir,
    //                                            "amrexinfo.raw");



    // std::vector<std::vector<double>> all_rmses(num_components);

    // iterator.iterate([&](int t, int l, int b) {
    //     const auto& actual = boxes[t][l][b];
    //     const auto& regen  = regen_boxes[t][l][b];

    //     std::vector<double> rmse = calc_rmse_per_box(actual, regen, num_components);

    //     for (int c = 0; c < num_components; ++c) {
    //         all_rmses[c].push_back(rmse[c]);
    //     }
    // });

    // // Compute mean RMSE and adjusted loss
    // for (int c = 0; c < num_components; ++c) {
    //     double mean_rmse = std::accumulate(all_rmses[c].begin(), all_rmses[c].end(), 0.0) /
    //                        all_rmses[c].size();

    //     spdlog::info("RMSE, {} = {}", amrexinfo_read.comp_names[c], mean_rmse);

    //     double loss = calc_adj_loss(mean_rmse, max_values[c] - min_values[c]);
    //     spdlog::info("Adjusted loss, {} = {}", amrexinfo_read.comp_names[c], loss);
    // }

    // LocDimData locs_read = read_loc_dim_from_bin(cfg.compressed_dir,
    //                                              "locations.raw",
    //                                              counts_read,
    //                                              iterator,
    //                                              num_times,
    //                                              num_levels);

    // LocDimData dims_read = read_loc_dim_from_bin(cfg.compressed_dir,
    //                                              "dimensions.raw",
    //                                              counts_read,
    //                                              iterator,
    //                                              num_times,
    //                                              num_levels);

    // write_plotfiles(regen_boxes,
    //                 locs_read,
    //                 dims_read,
    //                 num_times,
    //                 num_levels,
    //                 num_components,
    //                 amrexinfo_read,
    //                 "../../regenerated-plotfiles/");

    // spdlog::info("Sucessfully wrote plotfiles.");

//     amrex::Finalize();
// }

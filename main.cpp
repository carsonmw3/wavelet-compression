#include "preprocess.h"
#include "box-structs.h"
#include "readandwrite.h"
#include "compressor.h"
#include "decompressor.h"
#include "calc-loss.h"
#include "writeplotfile.h"

#include <spdlog/spdlog.h>
#include <chrono>

#include <AMReX.H>
#include <AMReX_ParmParse.H>



int main(int argc, char* argv[]) {

    // TODO: decide on input strategy
    amrex::Initialize(argc, argv);
    spdlog::set_level(spdlog::level::debug);

    std::string      data_dir;
    int              min_time;
    int              max_time;
    int              min_level;
    int              max_level;
    std::vector<int> components;
    float            keep;
    std::string      compressed_dir;
    amrex::ParmParse pp;

    pp.query("datadir", data_dir);
    pp.query("mintime", min_time);
    pp.query("maxtime", max_time);
    pp.query("minlevel", min_level);
    pp.query("maxlevel", max_level);
    pp.queryarr("components", components);
    pp.query("keep", keep);
    pp.query("compressedDir", compressed_dir);

    std::vector<std::string> files;
    for (int t = min_time; t < max_time+1; t++) {
        std::string filename = data_dir +"plt0" + std::to_string(t) + "00";
        files.push_back(filename);
    }

    std::vector<int> levels;
    for (int l = min_level; l < max_level+1; l++) {
        levels.push_back(l);
    }

    spdlog::info("Processing data...");
    auto start = std::chrono::high_resolution_clock::now();

    int num_times      = max_time - min_time + 1;
    int num_levels     = max_level - min_level + 1;
    int num_components = components.size();

    AllData data = preprocess_data(files,
                                   components,
                                   levels);

    auto& boxes       = data.boxes;
    auto& locations   = data.locations;
    auto& dimensions  = data.dimensions;
    auto& box_counts  = data.box_counts;
    auto& min_values  = data.min_values;
    auto& max_values  = data.max_values;
    auto& amrexinfo   = data.amrexinfo;

    write_loc_dim_to_bin(locations,
                         compressed_dir,
                         "locations.raw",
                         num_times,
                         num_levels);
    write_loc_dim_to_bin(dimensions,
                         compressed_dir,
                         "dimensions.raw",
                         num_times,
                         num_levels);
    write_box_counts(box_counts,
                     compressed_dir,
                     "boxcounts.raw",
                     num_times,
                     num_levels);

    auto end      = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end - start).count();
    spdlog::info(
        "Successfully processed data in {} seconds. Beginning compression...",
        duration);

    auto start1 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t <= max_time - min_time; t++) {
        for (int lev = 0; lev <= max_level - min_level; lev++) {
            for (int box_idx = 0; box_idx < boxes[t][lev].size(); box_idx++) {

                multiBox3D& current_box = boxes[t][lev][box_idx];
                compress(
                    current_box, components, keep, t, lev, box_idx, compressed_dir);
            }
        }
    }

    auto end1      = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration<double>(end1 - start1).count();
    spdlog::info(
        "Compression completed in {} seconds. Beginning decompression...",
        duration1);

    auto start2      = std::chrono::high_resolution_clock::now();
    auto counts_read = read_box_counts(compressed_dir,
                                       "boxcounts.raw",
                                       num_times,
                                       num_levels);

    std::vector<std::vector<std::vector<multiBox3D>>> regen_boxes;

    for (int t = 0; t <= max_time - min_time; t++) {
        std::vector<std::vector<multiBox3D>> regen_boxes_t;
        for (int lev = 0; lev <= max_level - min_level; lev++) {
            std::vector<multiBox3D> regen_boxes_l;
            for (int box_idx = 0; box_idx < boxes[t][lev].size(); box_idx++) {

                multiBox3D multibox;

                for (int component : components) {

                    std::string file_path = compressed_dir + "compressed-wavelet-" +
                                            std::to_string(t) + "-" +
                                            std::to_string(lev) + "-" +
                                            std::to_string(component) + "-" +
                                            std::to_string(box_idx) + ".xz";

                    Box3D curr_box = decompress(file_path, t, lev, component, box_idx);

                    multibox.push_back(std::move(curr_box));

                }
                regen_boxes_l.push_back(std::move(multibox));
            }
            regen_boxes_t.push_back(std::move(regen_boxes_l));
        }
        regen_boxes.push_back(std::move(regen_boxes_t));
    }

    auto end2      = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration<double>(end2 - start2).count();
    spdlog::info("Decompression completed in {} seconds. Calculating loss...",
                 duration2);


    std::vector<double> rmses = calc_avg_rmse(boxes, regen_boxes, num_components);
    for (int c = 0; c < num_components; c++) {
        spdlog::info("RMSE, {} = {}", amrexinfo.comp_names[c], rmses[c]);
        double loss = calc_adj_loss(rmses[c], max_values[c]-min_values[c]);
        spdlog::info("Adjusted loss, {} = {}", amrexinfo.comp_names[c], loss);
    }

    LocDimData locs_read = read_loc_dim_from_bin(compressed_dir,
                                                 "locations.raw",
                                                 counts_read,
                                                 num_times,
                                                 num_levels);

    LocDimData dims_read = read_loc_dim_from_bin(compressed_dir,
                                                 "dimensions.raw",
                                                 counts_read,
                                                 num_times,
                                                 num_levels);

    write_plotfiles(regen_boxes,
                    locs_read,
                    dims_read,
                    num_times,
                    num_levels,
                    num_components,
                    amrexinfo,
                    "../../regenerated-plotfiles/");

    spdlog::info("Sucessfully wrote plotfiles.");

    amrex::Finalize();
}

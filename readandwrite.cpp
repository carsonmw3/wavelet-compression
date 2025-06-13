#include "readandwrite.h"

#include <spdlog/spdlog.h>
#include <fstream>


// TODO: Note the addition of 'static'. This is a weird c++ thing; if you
// only use a function inside a single .cpp, it is recommended to tell the
// compiler that.
static void write_float(std::ofstream& stream, float val) {
    stream.write(reinterpret_cast<const char*>(&val), sizeof(float));
}


static float read_float(std::ifstream& stream) {
    float value;
    stream.read(reinterpret_cast<char*>(&value), sizeof(float));
    return value;
}


static std::ofstream open_write(std::string path, std::string outfile) {
    std::string   filename = path + outfile;
    std::ofstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        spdlog::error("Failed to open file: {}", filename);
        exit(EXIT_FAILURE);
    }

    return file;
}



static std::ifstream open_read(std::string path, std::string infile) {
    std::string   filename = path + infile;
    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        spdlog::error("Failed to open file: {}", filename);
        exit(EXIT_FAILURE);
    }

    return file;
}



void write_loc_dim_to_bin(LocDimData  data,
                          std::string path,
                          std::string out_file,
                          int         num_times,
                          int         num_levels) {

    std::ofstream file = open_write(path, out_file);

    for (int t = 0; t < num_times; t++) {

        for (int lev = 0; lev < num_levels; lev++) {

            for (int box_idx = 0; box_idx < data[t][lev].size(); box_idx++) {

                for (int coord = 0; coord < 3; coord++) {

                    float value = data[t][lev][box_idx][coord];
                    write_float(file, value);
                }
            }
        }
    }

    file.close();
}



LocDimData read_loc_dim_from_bin(std::string const&            path,
                                 std::string const&            in_file,
                                 std::vector<std::vector<int>> counts,
                                 int                           num_times,
                                 int                           num_levels) {

    std::ifstream file = open_read(path, in_file);

    LocDimData out;

    for (int t = 0; t < num_times; t++) {
        std::vector<std::vector<std::vector<int>>> out_t;

        for (int lev = 0; lev < num_levels; lev++) {
            std::vector<std::vector<int>> out_l;
            if (t >= counts.size()) {
                throw std::runtime_error("Time index out of bounds in counts vector");
            }
            if (lev >= counts[t].size()) {
                throw std::runtime_error("Level index out of bounds in counts vector");
            }
            int                           curr_num_boxes = counts[t][lev];

            for (int box_idx = 0; box_idx < curr_num_boxes; box_idx++) {
                std::vector<int> current_coords;

                for (int coord = 0; coord < 3; coord++) {
                    float value = read_float(file);
                    current_coords.push_back(value);
                }

                out_l.push_back(current_coords);
            }

            out_t.push_back(out_l);
        }

        out.push_back(out_t);
    }

    file.close();
    return out;
}


void write_box_counts(std::vector<std::vector<int>> counts,
                      std::string const&            path,
                      std::string const&            out_file,
                      int                           num_times,
                      int                           num_levels) {

    std::ofstream file = open_write(path, out_file);

    for (int t = 0; t < num_times; t++) {

        for (int lev = 0; lev < num_levels; lev++) {

            int value = counts.at(t).at(lev);
            write_float(file, value);
        }
    }

    file.close();
}


std::vector<std::vector<int>> read_box_counts(std::string path,
                                              std::string in_file,
                                              int         num_times,
                                              int         num_levels) {

    std::ifstream file = open_read(path, in_file);

    std::vector<std::vector<int>> read_counts;

    for (int t = 0; t < num_times; t++) {
        std::vector<int> read_counts_t;

        for (int lev = 0; lev < num_levels; lev++) {
            int value = static_cast<int>(read_float(file));
            read_counts_t.push_back(value);
        }

        read_counts.push_back(read_counts_t);
    }

    file.close();
    return read_counts;
}


void write_box_with_loc_dim(Volume3D           volume,
                            std::string        path,
                            std::string        out_file,
                            std::vector<float> location,
                            std::vector<float> dimension) {

    std::ofstream file = open_write(path, out_file);

    // Write location
    for (float coord : location) {
        write_float(file, coord);
    }

    // Write dimension
    for (float dim : dimension) {
        write_float(file, dim);
    }


    // Write volume data
    volume.iterate([&](float value, size_t, size_t, size_t) {
        write_float(file, value);
    });
    // TODO: Check and see if this can be used anywhere else

    file.close();
}


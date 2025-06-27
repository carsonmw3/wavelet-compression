#include "readandwrite.h"

#include <spdlog/spdlog.h>
#include <fstream>
#include <filesystem>

#include <doctest/doctest.h>

static void write_float(std::ofstream& stream, float val) {
    stream.write(reinterpret_cast<const char*>(&val), sizeof(float));
}


static float read_float(std::ifstream& stream) {
    float value;
    stream.read(reinterpret_cast<char*>(&value), sizeof(float));
    return value;
}

static void write_string(std::ofstream& stream, std::string& str) {
    size_t len = str.size();
    stream.write(reinterpret_cast<const char*>(&len), sizeof(len));
    stream.write(str.c_str(), len);
}

static std::string read_string(std::ifstream& stream) {
    size_t len;
    stream.read(reinterpret_cast<char*>(&len), sizeof(len));
    std::string str(len, '\0');
    stream.read(&str[0], len);
    return str;
}

static void write_int(std::ofstream& stream, int val) {
    stream.write(reinterpret_cast<const char*>(&val), sizeof(int));
}

static int read_int(std::ifstream& stream) {
    int value;
    stream.read(reinterpret_cast<char*>(&value), sizeof(int));
    return value;
}

static void write_long_double(std::ofstream& stream, long double val) {
    stream.write(reinterpret_cast<const char*>(&val), sizeof(long double));
}

static long double read_long_double(std::ifstream& stream) {
    long double value;
    stream.read(reinterpret_cast<char*>(&value), sizeof(long double));
    return value;
}

static void write_double(std::ofstream& stream, double val) {
    stream.write(reinterpret_cast<const char*>(&val), sizeof(double));
}

static double read_double(std::ifstream& stream) {
    double value;
    stream.read(reinterpret_cast<char*>(&value), sizeof(double));
    return value;
}

static void write_vector_string(std::ofstream& stream, const std::vector<std::string>& vec) {
    size_t size = vec.size();
    stream.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (auto& str : vec) {
        write_string(stream, const_cast<std::string&>(str));
    }
}

static std::vector<std::string> read_vector_string(std::ifstream& stream) {
    size_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(size));
    std::vector<std::string> vec(size);
    for (size_t i = 0; i < size; ++i) {
        vec[i] = read_string(stream);
    }
    return vec;
}

static void write_vector_int(std::ofstream& stream, const std::vector<int>& vec) {
    size_t size = vec.size();
    stream.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (int val : vec) {
        write_int(stream, val);
    }
}

static std::vector<int> read_vector_int(std::ifstream& stream) {
    size_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(size));
    std::vector<int> vec(size);
    for (size_t i = 0; i < size; ++i) {
        vec[i] = read_int(stream);
    }
    return vec;
}

static void write_vector_long_double(std::ofstream& stream, const std::vector<long double>& vec) {
    size_t size = vec.size();
    stream.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (long double val : vec) {
        write_long_double(stream, val);
    }
}

static std::vector<long double> read_vector_long_double(std::ifstream& stream) {
    size_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(size));
    std::vector<long double> vec(size);
    for (size_t i = 0; i < size; ++i) {
        vec[i] = read_long_double(stream);
    }
    return vec;
}

static void write_vector_vector_double(std::ofstream& stream, const std::vector<std::vector<double>>& vecvec) {
    size_t outer = vecvec.size();
    stream.write(reinterpret_cast<const char*>(&outer), sizeof(outer));
    for (const auto& vec : vecvec) {
        size_t inner = vec.size();
        stream.write(reinterpret_cast<const char*>(&inner), sizeof(inner));
        for (double val : vec) {
            write_double(stream, val);
        }
    }
}

static std::vector<std::vector<double>> read_vector_vector_double(std::ifstream& stream) {
    size_t outer;
    stream.read(reinterpret_cast<char*>(&outer), sizeof(outer));
    std::vector<std::vector<double>> vecvec(outer);
    for (size_t i = 0; i < outer; ++i) {
        size_t inner;
        stream.read(reinterpret_cast<char*>(&inner), sizeof(inner));
        std::vector<double> inner_vec(inner);
        for (size_t j = 0; j < inner; ++j) {
            inner_vec[j] = read_double(stream);
        }
        vecvec[i] = inner_vec;
    }
    return vecvec;
}

static void write_vector_vector_int(std::ofstream& stream, const std::vector<std::vector<int>>& vecvec) {
    size_t outer = vecvec.size();
    stream.write(reinterpret_cast<const char*>(&outer), sizeof(outer));
    for (const auto& vec : vecvec) {
        size_t inner = vec.size();
        stream.write(reinterpret_cast<const char*>(&inner), sizeof(inner));
        for (int val : vec) {
            write_int(stream, val);
        }
    }
}

static std::vector<std::vector<int>> read_vector_vector_int(std::ifstream& stream) {
    size_t outer;
    stream.read(reinterpret_cast<char*>(&outer), sizeof(outer));
    std::vector<std::vector<int>> vecvec(outer);
    for (size_t i = 0; i < outer; ++i) {
        size_t inner;
        stream.read(reinterpret_cast<char*>(&inner), sizeof(inner));
        std::vector<int> inner_vec(inner);
        for (size_t j = 0; j < inner; ++j) {
            inner_vec[j] = read_int(stream);
        }
        vecvec[i] = inner_vec;
    }
    return vecvec;
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
                          AMRIterator iterator) {

    std::ofstream file = open_write(path, out_file);

    iterator.iterate([&](int t, int lev, int box_idx) {
        for (int coord = 0; coord < 3; coord++) {
            float value = data[t][lev][box_idx][coord];
            write_float(file, value);
        }
    });

    file.close();
}



LocDimData read_loc_dim_from_bin(std::string const&            path,
                                 std::string const&            in_file,
                                 std::vector<std::vector<int>> counts,
                                 AMRIterator                   iterator,
                                                               int num_times,
                                                               int num_levels) {

    std::ifstream file = open_read(path, in_file);

    LocDimData out(num_times, std::vector<std::vector<std::vector<int>>>(num_levels));

    iterator.iterate([&](int t, int lev, int box_idx) {
        std::vector<int> current_coords;
        for (int coord = 0; coord < 3; ++coord) {
            float value = read_float(file);
            current_coords.push_back(value);
        }
        out[t][lev].push_back(std::move(current_coords));
    });

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


void write_amrexinfo(AMReXInfo   info,
                     std::string path,
                     std::string out_file) {

    std::ofstream file = open_write(path, out_file);

    write_vector_string(file, info.comp_names);
    write_vector_vector_double(file, info.geomcellinfo);
    write_vector_int(file, info.ref_ratios);
    write_vector_long_double(file, info.true_times);
    write_vector_vector_int(file, info.level_steps);

    write_int(file, info.xDim);
    write_int(file, info.yDim);
    write_int(file, info.zDim);

    file.close();

}


AMReXInfo read_amrex_info(std::string path,
                          std::string in_file) {

    std::ifstream file = open_read(path, in_file);
    AMReXInfo info;

    info.comp_names   = read_vector_string(file);
    info.geomcellinfo = read_vector_vector_double(file);
    info.ref_ratios   = read_vector_int(file);
    info.true_times   = read_vector_long_double(file);
    info.level_steps  = read_vector_vector_int(file);

    info.xDim = read_int(file);
    info.yDim = read_int(file);
    info.zDim = read_int(file);

    file.close();
    return info;
}


TEST_CASE("Read/write Loc/Dim data") {

    std::vector<int> testloc = { 0, 14, 44 };
    std::vector<std::vector<int>> vec;
    vec.push_back(testloc);
    std::vector<std::vector<std::vector<int>>> vec2;
    vec2.push_back(vec);
    vec2.push_back(vec);
    LocDimData test;
    test.push_back(vec2);
    test.push_back(vec2);

    std::string temp_dir = std::filesystem::temp_directory_path();

    AMRIterator iterator(2, 2, {{ 1, 1 }, { 1, 1 }}, 1);

    write_loc_dim_to_bin(test, temp_dir, "test.raw", iterator);
    LocDimData result = read_loc_dim_from_bin(temp_dir, "test.raw", { { 1, 1 }, { 1, 1 } }, iterator, 2, 2);

    REQUIRE(result == test);

}


TEST_CASE("Read/write Box counts") {

    std::vector<int> counts = { 403, 404, 333 };
    std::vector<std::vector<int>> test;
    test.push_back(counts);
    test.push_back(counts);

    std::string temp_dir = std::filesystem::temp_directory_path();

    write_box_counts(test, temp_dir, "test.raw", 2, 3);
    std::vector<std::vector<int>> result = read_box_counts(temp_dir, "test.raw", 2, 3);

    REQUIRE(result == test);

}


TEST_CASE("Read/write amrexinfo") {

    AMReXInfo test;

    test.comp_names = { "temp", "pressure" };
    test.geomcellinfo = { {0.6, 0.5, 0.4}, {0.8, 0.9, 1.0} };
    test.ref_ratios = { 2, 2, 2 };
    test.true_times = { 0.2219392, 0.3874982 };
    test.level_steps = { {1200, 1500}, {1800, 2000} };
    test.xDim = 256;
    test.yDim = 512;
    test.zDim = 256;

    std::string temp_dir = std::filesystem::temp_directory_path();

    write_amrexinfo(test, temp_dir, "test.raw");
    AMReXInfo result = read_amrex_info(temp_dir, "test.raw");

    REQUIRE(result.comp_names == test.comp_names);
    REQUIRE(result.geomcellinfo == test.geomcellinfo);
    REQUIRE(result.ref_ratios == test.ref_ratios);
    REQUIRE(result.true_times == test.true_times);
    REQUIRE(result.level_steps == test.level_steps);
    REQUIRE(result.xDim == test.xDim);
    REQUIRE(result.yDim == test.yDim);
    REQUIRE(result.zDim == test.zDim);

}




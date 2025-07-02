#include "argparse.h"

#include <AMReX_ParmParse.H>
#include <spdlog/spdlog.h>

#include <filesystem>

// parses user inputs in compression/estimate mode
Config parse_config_compress() {

    Config cfg;
    amrex::ParmParse pp;

    // Location of raw data (this dir contains "plt..." directories)
    pp.query("datadir", cfg.data_dir);

    // lowest timestep to compress (inclusive)
    pp.query("mintime", cfg.min_time);

    // highest timestep to compress (inclusive)
    pp.query("maxtime", cfg.max_time);

    // lowest level to compress (inclusive)
    pp.query("minlevel", cfg.min_level);

    // highest level to compress (inclusive)
    pp.query("maxlevel", cfg.max_level);

    // components to compress (e.g., for 6 and 25 type "components=6 25")
    pp.queryarr("components", cfg.components);

    // percent of wavelet coefficients that will be kept in compression.
    // higher keep value should lead to less compression but higher accuracy.
    // to start, try keep=0.99, 0.999, and 0.9999 in -estimate mode
    pp.query("keep", cfg.keep);

    // directory where the compressed data will be written (must already exist)
    pp.query("compresseddir", cfg.compressed_dir);

    return cfg;
}


// parses user inputs in decompression mode
Config parse_config_decompress() {
    Config cfg;
    amrex::ParmParse pp;

    // directory where the compressed data is stored
    pp.query("compresseddir", cfg.compressed_dir);

    return cfg;

}

// checks for a flag in the user input
bool has_flag(int argc, char* argv[], const std::string& flag) {
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == flag) {
            return true;
        }
    }
    return false;
}


// Removes all non-digit characters, as well as leading 0's from a
// filename
int clean_string(std::string filename) {

    std::string digits;

    // only keep digits
    for (char ch : filename) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            digits += ch;
        }
    }

    // get rid of leading 0's
    size_t start = digits.find_first_not_of('0');
    if (start == std::string::npos) {
        return 0;  // all zeros
    }
    digits = digits.substr(start);

    // convert to int
    return std::stoi(digits);

}


// infers vector of directories to be compressed based on inputs
std::vector<std::string> format_files(std::string data_dir,
                                      std::string min_time,
                                      std::string max_time) {

    int first = clean_string(min_time);
    int last = clean_string(max_time);

    std::vector<std::string> files;

    spdlog::info("This run involves the following files:");

    for (const auto& entry : std::filesystem::directory_iterator(data_dir)) {

        int current = clean_string(entry.path().string());

        if (current >= first && current <= last) {
            files.push_back(entry.path().string());
            spdlog::info("{}", entry.path().string());
        }

    }

    return files;

}


// constructs vector of levels to be compressed based on inputs
std::vector<int> format_levels(int min_level, int max_level) {

    std::vector<int> levels;
    for (int l = min_level; l <= max_level; ++l)
        levels.push_back(l);

    return levels;

}

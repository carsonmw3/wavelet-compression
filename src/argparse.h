#pragma once

#include <string>
#include <vector>

// stores all user inputs
struct Config {
    std::string              data_dir;
    std::string              compressed_dir;
    std::string              out_dir;
    std::string              min_time, max_time;
    int                      min_level, max_level;
    float                    keep;
    std::vector<std::string> components;
};

// parses user inputs in compression/estimate mode
Config parse_config_compress();

// parses user inputs in decompression mode
Config parse_config_decompress();

// checks for a flag in the user input
bool has_flag(int argc, char* argv[], const std::string& flag);

int clean_string(std::string filename);

// infers vector of directories to be compressed based on inputs
std::vector<std::string> format_files(std::string data_dir,
                                      std::string min_time,
                                      std::string max_time);

// constructs vector of levels to be compressed based on inputs
std::vector<int> format_levels(int min_level, int max_level);

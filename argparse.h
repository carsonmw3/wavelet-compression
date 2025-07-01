#pragma once

#include <string>
#include <vector>

// stores all user inputs
struct Config {
    std::string data_dir;
    std::string compressed_dir;
    int min_time, max_time;
    int min_level, max_level;
    float keep;
    std::vector<int> components;
};

// parses user inputs
Config parse_config();

#pragma once

#include "box-structs.h"

// Collect all relevant data for the specified compression run
AllData preprocess_data(std::vector<std::string> files,
                        std::vector<std::string> components,
                        std::vector<int>         levels);


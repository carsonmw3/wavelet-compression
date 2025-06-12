#pragma once

#include "box-structs.h"

AllData preprocess_data(std::vector<std::string> files,
                        std::vector<int>         components,
                        std::vector<int>         levels);

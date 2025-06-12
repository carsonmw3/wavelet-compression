#pragma once

#include "box-structs.h"

double calc_avg_rmse(std::vector<std::vector<std::vector<Box3D>>> const& actuals,
                     std::vector<std::vector<std::vector<Box3D>>> const& regens);

double calc_adj_loss(double rmse, double range);

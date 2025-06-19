#pragma once

#include "box-structs.h"

// TODO: We use this nested lookup structure quite a bit.
// Seems ripe for a class to encapsulate it. Seems like it could help out
// with iteration and processing later on
double calc_avg_rmse(std::vector<std::vector<std::vector<Box3D>>> const& actuals,
                     std::vector<std::vector<std::vector<Box3D>>> const& regens);

double calc_adj_loss(double rmse, double range);


#pragma once

#include "box-structs.h"

// TODO: We use this nested lookup structure quite a bit.
// Seems ripe for a class to encapsulate it. Seems like it could help out
// with iteration and processing later on
std::vector<double> calc_avg_rmse(std::vector<std::vector<std::vector<multiBox3D>>> const& actuals,
                                  std::vector<std::vector<std::vector<multiBox3D>>> const& regens,
                                  int                                                      num_components);

double calc_adj_loss(double rmse, double range);


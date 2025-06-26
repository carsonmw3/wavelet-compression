#pragma once

#include "box-structs.h"

std::vector<double> calc_avg_rmse(std::vector<std::vector<std::vector<multiBox3D>>> const& actuals,
                                  std::vector<std::vector<std::vector<multiBox3D>>> const& regens,
                                  int                                                      num_components);

double calc_adj_loss(double rmse, double range);

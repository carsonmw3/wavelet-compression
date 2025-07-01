#pragma once

#include "box-structs.h"

// Calculates root mean squared error for two multiboxes
std::vector<double> calc_rmse_per_box(const multiBox3D& actual,
                                      const multiBox3D& pred,
                                      int num_components);

// Component specific, calculates an adjusted loss metric that can
// be compared across components by dividing the rmse by the range
// of the data for that component
double calc_adj_loss(double rmse, double range);

// Calculates the size of the directory at path
double calc_size(std::string path);

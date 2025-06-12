#include "calc-loss.h"

#include <numeric>
#include <cmath>

double calc_avg_rmse(std::vector<std::vector<std::vector<Box3D>>> const& actuals,
              std::vector<std::vector<std::vector<Box3D>>> const& regens) {

    std::vector<double> rmses;

    for (int t = 0; t < actuals.size(); t++) {

        auto const& actuals_t = actuals[t];
        auto const& regens_t  = regens[t];

        for (int l = 0; l < actuals_t.size(); l++) {

            auto const& actuals_l = actuals_t[l];
            auto const& regens_l  = regens_t[l];

            for (int box_idx = 0; box_idx < actuals_l.size(); box_idx++) {

                auto const& actual = actuals_l[box_idx];
                auto const& pred   = regens_l[box_idx];

                int xdim = actual.width();
                int ydim = actual.height();
                int zdim = actual.depth();

                double sum = 0.0;

                for (int k = 0; k < zdim; k++) {
                    for (int j = 0; j < ydim; j++) {
                        for (int i = 0; i < xdim; i++) {
                            double diff = actual(i, j, k) - pred(i, j, k);
                            sum += diff * diff;
                        }
                    }
                }

                double rmse = sqrt(sum / (xdim * ydim * zdim));
                rmses.push_back(rmse);
            }
        }
    }

    double mean_rmse =
        accumulate(rmses.begin(), rmses.end(), 0.0) / rmses.size();
    return mean_rmse;
}


double calc_adj_loss(double rmse, double range) {
    return rmse / range;
}

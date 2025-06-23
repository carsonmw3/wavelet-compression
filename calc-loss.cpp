#include "calc-loss.h"

#include <numeric>
#include <cmath>
#include <spdlog/spdlog.h>

std::vector<double> calc_avg_rmse(std::vector<std::vector<std::vector<multiBox3D>>> const& actuals,
                     std::vector<std::vector<std::vector<multiBox3D>>> const& regens,
                     int                                                      num_components) {

    std::vector<std::vector<double>> rmses;
    rmses.resize(num_components);

    for (int t = 0; t < actuals.size(); t++) {

        auto const& actuals_t = actuals[t];
        auto const& regens_t  = regens[t];

        for (int l = 0; l < actuals_t.size(); l++) {

            auto const& actuals_l = actuals_t[l];
            auto const& regens_l  = regens_t[l];

            for (int box_idx = 0; box_idx < actuals_l.size(); box_idx++) {

                const auto& actual = actuals_l[box_idx];
                const auto& pred   = regens_l[box_idx];

                int xdim = actual[0].width();
                int ydim = actual[0].height();
                int zdim = actual[0].depth();

                std::vector<double> sums(num_components, 0.0);

                for (int c = 0; c < num_components; c++) {

                    const auto& actual_box = actual[c];
                    const auto& regen_box = pred[c];

                    if (!actual_box.is_valid()) {
                        spdlog::error("Invalid actual_box: dims = {}x{}x{}, expected = {}, actual = {}",
                                      actual_box.width(), actual_box.height(), actual_box.depth(),
                                      actual_box.width() * actual_box.height() * actual_box.depth(),
                                      actual_box.data_size());
                    }

                    for (int k = 0; k < zdim; k++) {
                        for (int j = 0; j < ydim; j++) {
                            for (int i = 0; i < xdim; i++) {
                                double diff = actual_box(i, j, k) - regen_box(i, j, k);
                                sums[c] += diff * diff;
                            }
                        }
                    }
                }

                for (int c = 0; c < num_components; c++) {
                    double rmse = sqrt(sums[c] / (xdim * ydim * zdim));
                    rmses[c].push_back(rmse);
                }
            }
        }
    }

    std::vector<double> mean_rmses;
    for (int c = 0; c < num_components; c++) {
        double mean_rmse =
            accumulate(rmses[c].begin(), rmses[c].end(), 0.0) / rmses[c].size();
        mean_rmses.push_back(mean_rmse);
    }
    return mean_rmses;
}


double calc_adj_loss(double rmse, double range) {
    return rmse / range;
}


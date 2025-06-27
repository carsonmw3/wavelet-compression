#include "calc-loss.h"

#include <cmath>
#include <spdlog/spdlog.h>

std::vector<double> calc_rmse_per_box(const multiBox3D& actual,
                                      const multiBox3D& pred,
                                      int num_components) {
    std::vector<double> rmse(num_components, 0.0);

    int xdim = actual[0].width();
    int ydim = actual[0].height();
    int zdim = actual[0].depth();

    for (int c = 0; c < num_components; c++) {
        const auto& actual_box = actual[c];
        const auto& regen_box  = pred[c];

        double sum = 0.0;

        for (int k = 0; k < zdim; ++k) {
            for (int j = 0; j < ydim; ++j) {
                for (int i = 0; i < xdim; ++i) {
                    double diff = actual_box(i, j, k) - regen_box(i, j, k);
                    sum += diff * diff;
                }
            }
        }

        rmse[c] = std::sqrt(sum / (xdim * ydim * zdim));
    }

    return rmse;
}


double calc_adj_loss(double rmse, double range) {
    return rmse / range;
}


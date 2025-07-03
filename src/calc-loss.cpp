#include "calc-loss.h"

#include <cmath>
#include <filesystem>
#include <spdlog/spdlog.h>

// Calculates root mean squared error for two multiboxes
// note: originally, iteration thru multiboxes happened within this function,
// but after implementing iterator.h, switched to this function calculating
// box by box
std::vector<double> calc_rmse_per_box(const multiBox3D& actual,
                                      const multiBox3D& pred,
                                      int num_components) {

    // storage for calculated rmses
    std::vector<double> rmse(num_components, 0.0);

    // extract box dimensions
    int xdim = actual[0].width();
    int ydim = actual[0].height();
    int zdim = actual[0].depth();

    for (int c = 0; c < num_components; c++) { // iterate thru components
        const auto& actual_box = actual[c];
        const auto& regen_box  = pred[c];

        double sum = 0.0; // total squared error

        for (int k = 0; k < zdim; ++k) {
            for (int j = 0; j < ydim; ++j) {
                for (int i = 0; i < xdim; ++i) {
                    double diff = actual_box(i, j, k) - regen_box(i, j, k); // error
                    sum += diff * diff; // squared error
                }
            }
        }

        rmse[c] = std::sqrt(sum / (xdim * ydim * zdim)); // root mean squared error
    }

    return rmse;
}


// Component specific, calculates an adjusted loss metric that can
// be compared across components by dividing the rmse by the range
// of the data for that component
double calc_adj_loss(double rmse, double range) {
    return rmse / range;
}


// Calculates the size of the directory at path
double calc_size(std::string path) {

    double size = 0;

    for (const auto& entry : std::filesystem::recursive_directory_iterator(path)) {
        size += std::filesystem::file_size(entry);
    }

    return size;

}


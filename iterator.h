#pragma once

#include <vector>

class AMRIterator {

private:
    size_t                               num_times;
    size_t                               num_levels;
    std::vector<std::vector<int>> const& box_counts;
    size_t                               num_components;

public:
    // Constructor
    AMRIterator(size_t                               num_times,
                size_t                               num_levels,
                std::vector<std::vector<int>> const& box_counts,
                size_t                               num_components)
        : num_times(num_times),
          num_levels(num_levels),
          box_counts(box_counts),
          num_components(num_components) { }

    template <class Function>
    void iterate(Function function) {
        for (int t = 0; t < num_times; t++) {
            for (int l = 0; l < num_levels; l++) {
                for (int b = 0; b < box_counts[t][l]; b++) {
                    function(t, l, b);
                }
            }
        }
    }
};

#include "preprocess.h"

#include <spdlog/spdlog.h>

#include <string>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_VisMF.H>
#include <AMReX_MFIter.H>


// Reads input data and stores it in a vector where the first dimension is the
// number of chunks of data (boxes), and each entry is a vector of floats where
// the first three are the location of the box, the next three are the dimensions of
// the box, and the rest is the data itself in order x, y, z.
std::tuple<std::vector<Box3D>,
           std::vector<Location>,
           std::vector<Dimensions>,
           int,
           float,
           float>
collectDataNewFormat (std::string lev_file,
                      int         component) {

    std::vector<Box3D>      boxes;
    std::vector<Location>   locations;
    std::vector<Dimensions> dimensions;

    int box_count = 0;

    float min_value = std::numeric_limits<float>::max();
    float max_value = std::numeric_limits<float>::min();

    amrex::MultiFab mf;

    // read data into multifab
    amrex::VisMF::Read(mf, lev_file);
    amrex::BoxArray ba = mf.boxArray();

    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.validbox();
        const amrex::Array4<amrex::Real>& mfdata = mf.array(mfi);

        const auto lo = lbound(box);
        const auto hi = ubound(box);

        const auto shape = box.size();

        Location loc;
        loc.push_back(lo.x);
        loc.push_back(lo.y);
        loc.push_back(lo.z);
        locations.push_back(loc);

        Dimensions dims;
        dims.push_back(shape[0]);
        dims.push_back(shape[1]);
        dims.push_back(shape[2]);
        dimensions.push_back(dims);

        Box3D box_data(shape[0], shape[1], shape[2]);

        for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    float value = mfdata(i,j,k,component);

                    box_data.set(i-lo.x, j-lo.y, k-lo.z, value);

                    if (value < min_value) {
                        min_value = value;
                    }

                    if (value > max_value) {
                        max_value = value;
                    }

                }
            }
        }

        boxes.push_back(std::move(box_data));
        box_count++;

    }
    
    mf.clear();

    return std::tuple(std::move(boxes),
                      locations,
                      dimensions,
                      box_count,
                      min_value,
                      max_value);
}


// TODO: add multi-component support
AllData preprocess_data(std::vector<std::string> files,
                           std::vector<int>         components,
                           std::vector<int>         levels) {

    AllData ret;

    auto& boxes      = ret.boxes;
    auto& locations  = ret.locations;
    auto& dimensions = ret.dimensions;
    auto& box_counts = ret.box_counts;
    auto& minval     = ret.min_value;
    auto& maxval     = ret.max_value;

    minval = std::numeric_limits<float>::max();
    maxval = std::numeric_limits<float>::min();

    for (int i = 0; i < files.size(); i++) {

        std::string filename = files[i];

        std::ifstream x;

        // open header
        std::string header = filename + "/Header";
        x.open(header.c_str(), std::ios::in);
        if (!x.is_open()) {
            spdlog::error("Failed to open header file: {}", filename);
        }

        // read in first line of header
        std::string str;
        x >> str;

        // read in number of components from header
        int nComp;
        x >> nComp;

        // read in variable names from header
        for (int n=0; n<nComp; n++) {
            x >> str;
        }

        // read in dimensionality from header
        int dim;
        x >> dim;

        if (dim != AMREX_SPACEDIM) {
            spdlog::error("Error: you are using a {}D build to open a {}D plotfile",
                          AMREX_SPACEDIM, dim);
        }

        std::vector<std::vector<Box3D>>      file_boxes;
        std::vector<std::vector<Location>>   file_locations;
        std::vector<std::vector<Dimensions>> file_dimensions;
        std::vector<int>                     file_box_counts;

        for (int level : levels) {

            std::string levX     = "/Level_"+std::to_string(level)+"/Cell";
            std::string lev_file = filename + levX;

            // TODO: add multi-component support
            std::tuple<std::vector<Box3D>,
                       std::vector<Location>,
                       std::vector<Dimensions>,
                       int,
                       float,
                       float>
                read_boxes = collectDataNewFormat(lev_file, components[0]);

	    spdlog::info("Processed data from time {}, level {}", i, level);
            file_boxes.push_back(std::move(std::get<0>(read_boxes)));
            file_locations.push_back(std::get<1>(read_boxes));
            file_dimensions.push_back(std::get<2>(read_boxes));
            file_box_counts.push_back(std::get<3>(read_boxes));

            float min = std::get<4>(read_boxes);
            if (min < minval) {
                minval = min;
            }

            float max = std::get<5>(read_boxes);
            if (max > maxval) {
                maxval = max;
            }

        }

        boxes.push_back(std::move(file_boxes));
        locations.push_back(file_locations);
        dimensions.push_back(file_dimensions);
        box_counts.push_back(file_box_counts);

    }

    return ret;

}


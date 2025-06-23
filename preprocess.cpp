#include "preprocess.h"

#include <spdlog/spdlog.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_VisMF.H>
#include <AMReX_MFIter.H>


// Reads input data and stores it in a vector where the first dimension is the
// number of chunks of data (boxes), and each entry is a vector of floats where
// the first three are the location of the box, the next three are the dimensions of
// the box, and the rest is the data itself in order x, y, z.
static LevelData collectDataNewFormat (std::string      lev_file,
                                       std::vector<int> components) {

    LevelData ret;

    auto& boxes      = ret.boxes;
    auto& locations  = ret.locations;
    auto& dimensions = ret.dimensions;
    auto& box_count  = ret.box_count;
    auto& min_values = ret.min_values;
    auto& max_values = ret.max_values;

    box_count = 0;
    min_values = std::vector<float>(components.size(),
                                    std::numeric_limits<float>::max());
    max_values = std::vector<float>(components.size(),
                                    std::numeric_limits<float>::min());

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

        multiBox3D box_data;

        for (int c = 0; c < components.size(); c++) {

            box_data.push_back(Box3D(shape[0], shape[1], shape[2], 0.0f));

            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {

                        float value = mfdata(i,j,k,components[c]);
                        box_data[c].set(i-lo.x, j-lo.y, k-lo.z, value);

                        if (value < min_values[c]) {
                            min_values[c] = value;
                        }

                        if (value > max_values[c]) {
                                max_values[c] = value;
                        }
                    }
                }
            }
        }

        boxes.push_back(std::move(box_data));
        box_count++;

    }

    mf.clear();

    return ret;
}


// TODO: add multi-component support
AllData preprocess_data(std::vector<std::string> files,
                        std::vector<int>         components,
                        std::vector<int>         levels) {

    AllData ret;

    auto& boxes        = ret.boxes;
    auto& locations    = ret.locations;
    auto& dimensions   = ret.dimensions;
    auto& box_counts   = ret.box_counts;
    auto& minvals      = ret.min_values;
    auto& maxvals      = ret.max_values;
    auto& amrexinfo    = ret.amrexinfo;

    minvals = std::vector<float>(components.size(),
                                    std::numeric_limits<float>::max());
    maxvals = std::vector<float>(components.size(),
                                    std::numeric_limits<float>::min());

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
        if (i == 0) {
            for (int n=0; n<nComp; n++) {
                x >> str;
                if (std::find(components.begin(),
                              components.end(),
                              n) != components.end()) {
                    amrexinfo.comp_names.push_back(str);
                }
            }
        } else {
            for (int n=0; n<nComp; n++) {
                x >> str;
            }
        }

        // read in dimensionality from header
        int dim;
        x >> dim;

        if (dim != AMREX_SPACEDIM) {
            spdlog::error("Error: you are using a {}D build to open a {}D plotfile",
                          AMREX_SPACEDIM, dim);
        }

        std::getline(x, str); // Consume remainder of the line after reading 'dim'

        // read in true time
        long double true_time;
        x >> true_time;
        amrexinfo.true_times.push_back(true_time);
        std::getline(x, str); // rest of line

        std::getline(x, str); // skip no. of levels

        // read in physical domain info
        std::vector<double> geomcell(6);
        std::getline(x, str);
        std::istringstream iss(str);
        double val1, val2, val3;
        iss >> val1 >> val2 >> val3;
        geomcell[0] = val1;
        geomcell[1] = val2;
        geomcell[2] = val3;

        std::getline(x, str);
        std::istringstream iss1(str);
        double val4, val5, val6;
        iss1 >> val4 >> val5 >> val6;
        geomcell[3] = val4;
        geomcell[4] = val5;
        geomcell[5] = val6;

        amrexinfo.geomcellinfo.push_back(geomcell);

        // read in refinement ratio
        if (i == 0) {
            std::vector<int> refratio(dim);
            std::getline(x, str);
            std::istringstream iss2(str);
            for (int i=0; i < dim; i++) {
                int ref;
                iss2 >> ref;
                refratio[i] = ref;
            }
            amrexinfo.ref_ratios = refratio;
        } else {
            std::getline(x, str);
        }

        // read in level dimensions
        std::getline(x, str);
        size_t first = str.find('(');
        first = str.find('(', first + 1);
        first = str.find('(', first + 1); // bounds we want are after third "("
        size_t end = str.find(')', first);

        std::string dims_str = str.substr(first+1, end-first+1);

        // split, convert to ints
        std::istringstream iss3(dims_str);
        std::string val;
        std::vector<int> dims;

        while (std::getline(iss3, val, ',')) {
            dims.push_back(std::stoi(val));
        }

        amrexinfo.xDim = dims[0] + 1;
        amrexinfo.yDim = dims[1] + 1;
        amrexinfo.zDim = dims[2] + 1;


        // read in level_steps
        std::vector<int> level_steps_i(levels.size());
        std::getline(x, str);
        std::istringstream iss4(str);
        for (int i=0; i < levels.size(); i++) {
            int ls;
            iss4 >> ls;
            level_steps_i[i] = ls;
        }
        amrexinfo.level_steps.push_back(level_steps_i);

        std::vector<std::vector<multiBox3D>> file_boxes;
        std::vector<std::vector<Location>>   file_locations;
        std::vector<std::vector<Dimensions>> file_dimensions;
        std::vector<int>                     file_box_counts;

        for (int level : levels) {

            std::string levX     = "/Level_"+std::to_string(level)+"/Cell";
            std::string lev_file = filename + levX;

            // TODO: add multi-component support
            LevelData read_boxes = collectDataNewFormat(lev_file, components);

            spdlog::info("Processed data from time {}, level {}", i, level);

            file_boxes.push_back(std::move(read_boxes.boxes));
            file_locations.push_back(read_boxes.locations);
            file_dimensions.push_back(read_boxes.dimensions);
            file_box_counts.push_back(read_boxes.box_count);

            for (int c = 0; c < components.size(); c++){

                float min = read_boxes.min_values[c];
                if (min < minvals[c]) {
                    minvals[c] = min;
                }

                float max = read_boxes.max_values[c];
                if (max > maxvals[c]) {
                    maxvals[c] = max;
                }
            }

        }

        boxes.push_back(std::move(file_boxes));
        locations.push_back(file_locations);
        dimensions.push_back(file_dimensions);
        box_counts.push_back(file_box_counts);

    }

    return ret;

}





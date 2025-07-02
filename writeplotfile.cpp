#include "writeplotfile.h"
#include "tmpdir.h"

#include <spdlog/spdlog.h>
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <AMReX_MultiFab.H>
#include <AMReX_BoxList.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>


// initializes a MultiFab with the proper layout of boxes according to inputs
static amrex::MultiFab initializeMF (std::vector<Location>   locations,
                                     std::vector<Dimensions> dimensions,
                                     int                     num_components) {

    amrex::BoxList boxes;

    // for every location, create the corresponding box
    for (int i = 0; i < locations.size(); i++) {

        Location   loc  = locations[i];
        Dimensions dims = dimensions[i];

        amrex::IntVect lo(loc[0], loc[1], loc[2]);
        amrex::IntVect hi(loc[0] + dims[0] - 1,
                          loc[1] + dims[1] - 1,
                          loc[2] + dims[2] - 1);

        amrex::Box current(lo, hi);

        boxes.push_back(current);

    }

    // create box array, distribution mapping from box list
    amrex::BoxArray            array(boxes);
    amrex::DistributionMapping dm(array);

    // no ghost cells
    int ngrow = 0;

    // create multifab
    amrex::MultiFab output(array, dm, num_components, ngrow);

    if (output.size() == 0) {
        spdlog::error("Error: MultiFab not initialized correctly.");
    }

    return output;

}


// take an existing multifab and populate it with data
static void populateMF (amrex::MultiFab&         multi,
                        std::vector<multiBox3D>& data,
                        int                      num_components) {

    // keep track of what box we're on for indexing into data
    int box_idx = 0;

    for (amrex::MFIter mfi(multi, false); mfi.isValid(); ++mfi) { // iterate thru each box in multifab

        if (box_idx >= data.size()) {
            spdlog::error("Index out of bounds: box_idx = {}, data.size() = {}", box_idx, data.size());
            std::abort();
        }

        // extract the box info for data entry
        const amrex::Box&                box    = mfi.validbox();
        const amrex::Array4<amrex::Real> mfdata = multi.array(mfi);
        const auto                       lo     = lbound(box);
        const auto                       hi     = ubound(box);

        // get the appropriate multibox of data
        multiBox3D& current = data[box_idx];

        for (int c = 0; c < num_components; c++) { // iterate thru components

            // get the box for the current component
            Box3D& curr_box = current[c];

            // populate multifab box with data
            for (int k = lo.z; k <= hi.z; k++) {

                for (int j = lo.y; j <= hi.y; j++) {

                    for (int i = lo.x; i <= hi.x; i++) {

                            mfdata(i,j,k,c) = curr_box.get(i-lo.x, j-lo.y, k-lo.z);

                    }
                }
            }
        }

        box_idx++;

    }

}


// writes a plotfile for each timestep of a compression run to the directory "out"
void write_plotfiles(std::vector<std::vector<std::vector<multiBox3D>>> &data,
                     LocDimData                                        locations,
                     LocDimData                                        dimensions,
                     int                                               num_times,
                     int                                               num_levels,
                     int                                               num_components,
                     std::vector<std::string>                          comp_names,
                     AMReXInfo                                         amrexinfo,
                     std::string                                       out) {

    for (int t = 0; t < num_times; t++) { // iterate thru timesteps

        // name of plotfile directory
        const std::string name = out + amrex::Concatenate("plt", t+74);

        // vectors to store plotfile info
        std::vector<amrex::MultiFab>   mfs;
        amrex::Vector<amrex::Geometry> geoms;
        amrex::Real                    time = amrexinfo.true_times[t];
        amrex::Vector<int>             level_steps_amr;
        amrex::Vector<amrex::IntVect>  ref_ratio;
        amrex::Vector<std::string>     varnames;

        // get names of compressed components
        for (std::string& comp : comp_names) {
            varnames.push_back(comp);
        }

        for (int l = 0; l < num_levels; l++) { // iterate thru levels

            std::vector<Location>    current_locs = locations[t][l];
            std::vector<Dimensions>  current_dims = dimensions[t][l];
            std::vector<multiBox3D>& current_data = data[t][l];

            // initialize and populate a multifab for this level
            amrex::MultiFab mf = initializeMF(current_locs, current_dims, num_components);
            populateMF(mf, current_data, num_components);

            mfs.push_back(std::move(mf));

            // increment the dimensions according to the refinement ratio
            int xDimCurrent = amrexinfo.xDim * pow(amrexinfo.ref_ratios[0], l);
            int yDimCurrent = amrexinfo.yDim * pow(amrexinfo.ref_ratios[1], l);
            int zDimCurrent = amrexinfo.zDim * pow(amrexinfo.ref_ratios[2], l);

            // specify the indexing domain
            amrex::Box domain(amrex::IntVect(0, 0, 0),
                              amrex::IntVect(xDimCurrent-1, yDimCurrent-1, zDimCurrent-1));

            // specify the physical domain
            std::vector<double> geomcell = amrexinfo.geomcellinfo[t];
            amrex::RealBox cell({geomcell[0], geomcell[1], geomcell[2]},
                                {geomcell[3], geomcell[4], geomcell[5]});

            // non-periodic
            amrex::Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0, 0, 0)};

            // create the geometry object for the plotfile (0 is for carteisan coodinate system)
            const amrex::Geometry geom(domain, cell, 0, is_periodic);
            geoms.push_back(geom);

            // time step associated with data
            level_steps_amr.push_back(amrexinfo.level_steps[t][l]);

            // refinement ratio
            if (l > 0) {
                amrex::IntVect ratio(amrexinfo.ref_ratios[0],
                                     amrexinfo.ref_ratios[1],
                                     amrexinfo.ref_ratios[2]);
                ref_ratio.push_back(ratio);
            }

        }

        // format vectors for writing function
        amrex::Vector<const amrex::MultiFab*> mfPtrs;

        for (auto& mf : mfs) {
            mfPtrs.push_back(&mf);
        }

        const amrex::Vector<const amrex::MultiFab*> constMfs        = mfPtrs;
        const amrex::Vector<std::string>            constVarnames   = varnames;
        const amrex::Vector<amrex::Geometry>        constGeoms      = geoms;
        const amrex::Vector<int>                    constLevelSteps = level_steps_amr;
        const amrex::Vector<amrex::IntVect>         constRefRatio   = ref_ratio;

        // write plotfile
        amrex::WriteMultiLevelPlotfile(name,
                                       num_levels,
                                       constMfs,
                                       constVarnames,
                                       constGeoms,
                                       time,
                                       constLevelSteps,
                                       constRefRatio);

    }

}


// for use in testing, checks if two files are identical
static bool files_are_identical(const std::string& file1, const std::string& file2) {
    std::ifstream f1(file1, std::ios::binary);
    std::ifstream f2(file2, std::ios::binary);

    if (!f1.is_open() || !f2.is_open()) return false;

    std::istreambuf_iterator<char> begin1(f1), end1;
    std::istreambuf_iterator<char> begin2(f2), end2;

    return std::equal(begin1, end1, begin2, end2);
}



TEST_CASE("Initialize, populate MultiFab") {

    int argc = 1;
    const char* arg0 = "testprog";
    char* argv0 = const_cast<char*>(arg0);
    char** argv = &argv0;

    amrex::Initialize(argc, argv);

    std::vector<Location>   test_locs      = { { 0, 0, 0 },
                                               { 16, 32, 64 } };
    std::vector<Dimensions> test_dims      = { { 16, 32, 64 },
                                               { 8, 4, 2 } };
    int                     num_components = 2;

    amrex::MultiFab test = initializeMF(test_locs, test_dims, num_components);

    REQUIRE(test.size() == test_locs.size());

    REQUIRE(test.isDefined());

    Box3D testbox1(16, 32, 64, 3902.4f);
    Box3D testbox2(8, 4, 2, 16.00f);

    multiBox3D test1;
    test1.push_back(testbox1.clone());
    test1.push_back(testbox1.clone());
    multiBox3D test2;
    test2.push_back(testbox2.clone());
    test2.push_back(testbox2.clone());
    std::vector<multiBox3D> testboxes;
    testboxes.push_back(std::move(test1));
    testboxes.push_back(std::move(test2));

    populateMF(test, testboxes, num_components);
    INFO("num boxes" << test.size());

    REQUIRE(test.size() == 2);

    amrex::Array4<amrex::Real> data1 = test.array(0);
    amrex::Array4<amrex::Real> data2 = test.array(1);

    const amrex::Box& box0 = test.boxArray()[0];
    const amrex::Box& box1 = test.boxArray()[1];

    auto lo0 = lbound(box0);
    auto lo1 = lbound(box1);

    float val1 = data1(lo0.x + 5, lo0.y + 10, lo0.z + 20, 0);
    float val2 = data2(lo1.x + 1, lo1.y + 1, lo1.z + 0, 1);

    REQUIRE(val1 == doctest::Approx(3902.4f).epsilon(0.001));
    REQUIRE(val2 == doctest::Approx(16.00f).epsilon(0.001));

    test.clear();

    amrex::Finalize();

}


TEST_CASE("Writing plotfiles") {


    int argc = 1;
    const char* arg0 = "testprog";
    char* argv0 = const_cast<char*>(arg0);
    char** argv = &argv0;

    amrex::Initialize(argc, argv);

    std::vector<Location>   test_locs      = { { 0, 0, 0 },
                                               { 16, 32, 64 } };
    std::vector<Dimensions> test_dims      = { { 16, 32, 64 },
                                               { 8, 4, 2 } };
    int                     num_components = 2;
    int                     num_times      = 2;
    int                     num_levels     = 2;

    Box3D testbox1(16, 32, 64, 3902.4f);
    Box3D testbox2(8, 4, 2, 16.00f);

    std::vector<std::vector<std::vector<multiBox3D>>> testdata;
    LocDimData locs;
    LocDimData dims;
    for (int i=0; i < num_times; i++) {

        std::vector<std::vector<multiBox3D>> test4;
        std::vector<std::vector<Location>> veclocs;
        std::vector<std::vector<Dimensions>> vecdims;
        for (int i=0; i < num_levels; i++) {

            multiBox3D test1;
            multiBox3D test2;

            for (int i=0; i < num_components; i++) {
                test1.push_back(testbox1.clone());
                test2.push_back(testbox2.clone());
            }

            std::vector<multiBox3D> test3;
            test3.push_back(std::move(test1));
            test3.push_back(std::move(test2));

            test4.push_back(std::move(test3));
            veclocs.push_back(test_locs);
            vecdims.push_back(test_dims);
        }

        testdata.push_back(std::move(test4));
        locs.push_back(veclocs);
        dims.push_back(vecdims);
    }

    AMReXInfo info;

    info.geomcellinfo = { {0.6, 0.5, 0.4, 0.8, 0.9, 1.0},
                          {0.6, 0.5, 0.4, 0.8, 0.9, 1.0} };
    info.ref_ratios = { 2, 2, 2 };
    info.true_times = { 0.2219392, 0.3874982 } ;
    info.level_steps = { {1200, 1500}, {1800, 2000} };
    info.xDim = 256;
    info.yDim = 512;
    info.zDim = 256;

    std::vector<std::string> comp_names = { "temp", "pressure" };

    TempDir scratch_dir;

    write_plotfiles(testdata,
                    locs,
                    dims,
                    num_times,
                    num_levels,
                    num_components,
                    comp_names,
                    info,
                    scratch_dir.path().string() + "/");

    REQUIRE(files_are_identical("../tests/plt00074/", scratch_dir.path() / "plt00074/"));

    amrex::Finalize();

}

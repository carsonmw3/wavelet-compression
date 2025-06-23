#include "writeplotfile.h"

#include <spdlog/spdlog.h>

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


static amrex::MultiFab initializeMF (std::vector<Location>   locations,
                                     std::vector<Dimensions> dimensions,
                                     int                     num_components) {

    amrex::BoxList boxes;

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

    amrex::BoxArray            array(boxes);
    amrex::DistributionMapping dm(array);

    int ngrow = 0;

    amrex::MultiFab output(array, dm, num_components, ngrow);

    if (output.size() == 0) {
        spdlog::error("Error: MultiFab not initialized correctly.");
    }

    return output;

}


static void populateMF (amrex::MultiFab&         multi,
                        std::vector<multiBox3D>& data,
                        int                      num_components) {

    int box_idx = 0;

    for (amrex::MFIter mfi(multi, false); mfi.isValid(); ++mfi) {

        if (box_idx >= data.size()) {
            spdlog::error("Index out of bounds: box_idx = {}, data.size() = {}", box_idx, data.size());
            std::abort();
        }

        const amrex::Box&                box    = mfi.validbox();
        const amrex::Array4<amrex::Real> mfdata = multi.array(mfi);
        const auto                       lo     = lbound(box);
        const auto                       hi     = ubound(box);

        multiBox3D& current = data[box_idx];

        for (int c = 0; c < num_components; c++) {

            Box3D& curr_box = current[c];

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


// TODO: generalize output parameters
void write_plotfiles(std::vector<std::vector<std::vector<multiBox3D>>> &data,
                     LocDimData                                   locations,
                     LocDimData                                   dimensions,
                     int                                          num_times,
                     int                                          num_levels,
                     int                                          num_components,
                     AMReXInfo                                    amrexinfo,
                     std::string                                  out) {

    for (int t = 0; t < num_times; t++) {

        const std::string name = out + amrex::Concatenate("plt", t+74);
        std::vector<amrex::MultiFab> mfs;
        amrex::Vector<amrex::Geometry> geoms;
        amrex::Real time = amrexinfo.true_times[t];
        amrex::Vector<int> level_steps_amr;
        amrex::Vector<amrex::IntVect> ref_ratio;
        amrex::Vector<std::string> varnames;

        for (std::string& comp : amrexinfo.comp_names) {
            varnames.push_back(comp);
        }

        for (int l = 0; l < num_levels; l++) {

            std::vector<Location>   current_locs = locations[t][l];
            std::vector<Dimensions> current_dims = dimensions[t][l];
            std::vector<multiBox3D>&     current_data = data[t][l];

            amrex::MultiFab mf = initializeMF(current_locs, current_dims, num_components);
            populateMF(mf, current_data, num_components);

            mfs.push_back(std::move(mf));

            int xDimCurrent = amrexinfo.xDim * pow(2, l);
            int yDimCurrent = amrexinfo.yDim * pow(2, l);
            int zDimCurrent = amrexinfo.zDim * pow(2, l);

            amrex::Box domain(amrex::IntVect(0, 0, 0),
                              amrex::IntVect(xDimCurrent-1, yDimCurrent-1, zDimCurrent-1));

            std::vector<double> geomcell = amrexinfo.geomcellinfo[t];
            amrex::RealBox cell({geomcell[0], geomcell[1], geomcell[2]},
                                {geomcell[3], geomcell[4], geomcell[5]});

            amrex::Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0, 0, 0)};

            const amrex::Geometry geom(domain, cell, 0, is_periodic);
            geoms.push_back(geom);

            level_steps_amr.push_back(amrexinfo.level_steps[t][l]);

            if (l > 0) {
                amrex::IntVect ratio(amrexinfo.ref_ratios[0],
                                     amrexinfo.ref_ratios[1],
                                     amrexinfo.ref_ratios[2]);
                ref_ratio.push_back(ratio);
            }

        }

        amrex::Vector<const amrex::MultiFab*> mfPtrs;

        for (auto& mf : mfs) {
            mfPtrs.push_back(&mf);
        }

        const amrex::Vector<const amrex::MultiFab*> constMfs        = mfPtrs;
        const amrex::Vector<std::string>            constVarnames   = varnames;
        const amrex::Vector<amrex::Geometry>        constGeoms      = geoms;
        const amrex::Vector<int>                    constLevelSteps = level_steps_amr;
        const amrex::Vector<amrex::IntVect>         constRefRatio   = ref_ratio;

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


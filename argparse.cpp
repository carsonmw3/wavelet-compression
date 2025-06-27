#include "argparse.h"

#include <AMReX_ParmParse.H>
#include <spdlog/spdlog.h>

Config parse_config() {
    Config cfg;
    amrex::ParmParse pp;
    pp.query("datadir", cfg.data_dir);
    pp.query("mintime", cfg.min_time);
    pp.query("maxtime", cfg.max_time);
    pp.query("minlevel", cfg.min_level);
    pp.query("maxlevel", cfg.max_level);
    pp.queryarr("components", cfg.components);
    pp.query("keep", cfg.keep);
    pp.query("compressedDir", cfg.compressed_dir);

    return cfg;
}

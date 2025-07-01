#include "argparse.h"

#include <AMReX_ParmParse.H>
#include <spdlog/spdlog.h>

// parses user inputs in compression/estimate mode
Config parse_config_compress() {
    Config cfg;
    amrex::ParmParse pp;

    // Location of raw data (this dir contains "plt..." directories)
    pp.query("datadir", cfg.data_dir);

    // lowest timestep to compress (inclusive)
    pp.query("mintime", cfg.min_time);

    // highest timestep to compress (inclusive)
    pp.query("maxtime", cfg.max_time);

    // lowest level to compress (inclusive)
    pp.query("minlevel", cfg.min_level);

    // highest level to compress (inclusive)
    pp.query("maxlevel", cfg.max_level);

    // components to compress (e.g., for 6 and 25 type "components=6 25")
    pp.queryarr("components", cfg.components);

    // percent of wavelet coefficients that will be kept in compression.
    // higher keep value should lead to less compression but higher accuracy.
    // to start, try keep=0.99, 0.999, and 0.9999 in -estimate mode
    pp.query("keep", cfg.keep);

    // directory where the compressed data will be written (must already exist)
    pp.query("compressedDir", cfg.compressed_dir);

    return cfg;
}


// parses user inputs in decompression mode
Config parse_config_decompress() {
    Config cfg;
    amrex::ParmParse pp;

    // directory where the compressed data is stored
    pp.query("compressedDir", cfg.compressed_dir);

    return cfg;

}

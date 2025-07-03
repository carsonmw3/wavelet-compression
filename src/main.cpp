#include "argparse.h"
#include "modes.h"

#include <spdlog/spdlog.h>

#include <AMReX.H>
#include <AMReX_ParmParse.H>


int main(int argc, char* argv[]) {

    amrex::Initialize(argc, argv);
    spdlog::set_level(spdlog::level::debug);

    if (has_flag(argc, argv, "-c")) {
        Config cfg = parse_config_compress();
        compress(cfg);
    } else if (has_flag(argc, argv, "-estimate")) {
        Config cfg = parse_config_compress();
        estimate(cfg);
    } else {
        Config cfg = parse_config_decompress();
        decompress(cfg);
    }

    amrex::Finalize();
    return 0;

}

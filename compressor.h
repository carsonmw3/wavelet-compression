#pragma once

#include "box-structs.h"

CompressedWavelet compress(Box3D const& box,
                           double       keep,
                           int          time,
                           int          level,
                           int          box_index,
                           std::string  compressed_dir);


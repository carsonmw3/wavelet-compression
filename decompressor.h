#pragma once

#include "box-structs.h"

Box3D decompress (std::string file_path,
                  int         time,
                  int         level,
                  int         component,
                  int         box_idx);


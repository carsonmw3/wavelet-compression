#include "grid.h"

#include <lzma.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

std::string compressed_dir = "../../wavelet/";

// Typedefs for readability
using Box3D    = Grid3D<float>;
using Volume3D = Grid3D<float>;
using LocDimData = std::vector<std::vector<std::vector<std::vector<int>>>>;

struct CompressedWavelet {
    std::vector<int> shape;
    std::vector<int> coeff_shape;
    std::vector<std::pair<int, int16_t>> rle_encoded;
};

struct DataResult {
    std::vector<std::vector<std::vector<Box3D>>> boxes;
    LocDimData                                   locations;
    LocDimData                                   dimensions;
    std::vector<std::vector<int>>                box_counts;
    float                                        min_value;
    float                                        max_value;
};


void write_float(std::ofstream& stream, float val) {
    stream.write(reinterpret_cast<const char*>(&val), sizeof(float));
}


float read_float(std::ifstream& stream) {
    float value;
    stream.read(reinterpret_cast<char*>(&value), sizeof(float));
    return value;
}



DataResult process_data(std::string const& dir_prefix,
                        int                min_time,
                        int                max_time,
                        int                min_level,
                        int                max_level,
                        int                component) {

    DataResult ret;

    auto& boxes      = ret.boxes;
    auto& locations  = ret.locations;
    auto& dimensions = ret.dimensions;
    auto& box_counts = ret.box_counts;
    auto& minval     = ret.min_value;
    auto& maxval     = ret.max_value;

    minval = std::numeric_limits<float>::max();
    maxval = std::numeric_limits<float>::min();

    for (int t = min_time; t <= max_time; ++t) {
        std::vector<std::vector<Box3D>>            boxes_t;
        std::vector<std::vector<std::vector<int>>> locations_t;
        std::vector<std::vector<std::vector<int>>> dimensions_t;
        std::vector<int>                           box_counts_t;

        for (int l = min_level; l <= max_level; ++l) {
            std::string filename = dir_prefix + "-" + std::to_string(l) + "/" +
                                   std::to_string(t) + "-wholeNewFormat-" +
                                   std::to_string(component) + "-" +
                                   std::to_string(l) + ".raw";

            std::ifstream file(filename, std::ios::binary);
            if (!file.is_open()) {
                spdlog::error("Could not open file: {}", filename);
                continue;
            }

            int                           box_count = 0;
            std::vector<Box3D>            boxes_l;
            std::vector<std::vector<int>> locations_l;
            std::vector<std::vector<int>> dimensions_l;

            while (file.peek() != EOF) {
                float x_f;
                file.read(reinterpret_cast<char*>(&x_f), sizeof(float));
                if (file.gcount() < sizeof(float)) break;
                float y_f = read_float(file);
                float z_f = read_float(file);

                std::vector<int> current_loc;
                int              x = static_cast<int>(x_f);
                current_loc.push_back(x);
                int y = static_cast<int>(y_f);
                current_loc.push_back(y);
                int z = static_cast<int>(z_f);
                current_loc.push_back(z);
                locations_l.push_back(current_loc);

                std::vector<int> current_dim;
                int              xdim = static_cast<int>(read_float(file));
                current_dim.push_back(x);
                int ydim = static_cast<int>(read_float(file));
                current_dim.push_back(y);
                int zdim = static_cast<int>(read_float(file));
                current_dim.push_back(z);
                dimensions_l.push_back(current_dim);

                Box3D box(xdim, ydim, zdim);

                for (int k = 0; k < zdim; ++k) {
                    for (int j = 0; j < ydim; ++j) {
                        for (int i = 0; i < xdim; ++i) {
                            float value  = read_float(file);
                            box(i, j, k) = value;

                            if (value < minval) {
                                minval = value;
                            }

                            if (value > maxval) {
                                maxval = value;
                            }
                        }
                    }
                }

                boxes_l.push_back(std::move(box));
                ++box_count;
            }

            boxes_t.push_back(std::move(boxes_l));
            locations_t.push_back(locations_l);
            dimensions_t.push_back(dimensions_l);
            box_counts_t.push_back(box_count);
        }

        boxes.push_back(std::move(boxes_t));
        locations.push_back(locations_t);
        dimensions.push_back(dimensions_t);
        box_counts.push_back(box_counts_t);
    }

    return ret;
}


std::ofstream open_write(std::string path, std::string outfile) {
    std::string   filename = path + outfile;
    std::ofstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        spdlog::error("Failed to open file: {}", filename);
        exit(EXIT_FAILURE);
    }

    return file;
}



std::ifstream open_read(std::string path, std::string infile) {
    std::string   filename = path + infile;
    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        spdlog::error("Failed to open file: {}", filename);
        exit(EXIT_FAILURE);
    }

    return file;
}



void write_loc_dim_to_bin(LocDimData  data,
                          std::string path,
                          std::string out_file,
                          int         max_time,
                          int         min_time,
                          int         min_level,
                          int         max_level) {

    std::ofstream file = open_write(path, out_file);

    for (int t = 0; t < max_time - min_time + 1; t++) {

        for (int lev = min_level; lev <= max_level; lev++) {

            for (int box_idx = 0; box_idx < data[t][lev].size(); box_idx++) {

                for (int coord = 0; coord < 3; coord++) {

                    float value = data[t][lev][box_idx][coord];
                    write_float(file, value);
                }
            }
        }
    }

    file.close();
}



LocDimData read_loc_dim_from_bin(std::string const&            path,
                                 std::string const&            in_file,
                                 std::vector<std::vector<int>> counts,
                                 int                           max_time,
                                 int                           min_time,
                                 int                           min_level,
                                 int                           max_level) {

    std::ifstream file = open_read(path, in_file);

    LocDimData out;

    for (int t = 0; t < max_time - min_time + 1; t++) {
        std::vector<std::vector<std::vector<int>>> out_t;

        for (int lev = min_level; lev <= max_level; lev++) {
            std::vector<std::vector<int>> out_l;
            int                           curr_num_boxes = counts[t][lev];

            for (int box_idx = 0; box_idx < curr_num_boxes; box_idx++) {
                std::vector<int> current_coords;

                for (int coord = 0; coord < 3; coord++) {
                    float value = read_float(file);
                    current_coords.push_back(value);
                }

                out_l.push_back(current_coords);
            }

            out_t.push_back(out_l);
        }

        out.push_back(out_t);
    }

    file.close();
    return out;
}


void write_box_counts(std::vector<std::vector<int>> counts,
                      std::string const&            path,
                      std::string const&            out_file,
                      int                           max_time,
                      int                           min_time,
                      int                           min_level,
                      int                           max_level) {

    std::ofstream file = open_write(path, out_file);

    for (int t = 0; t < max_time - min_time + 1; t++) {

        for (int lev = min_level; lev <= max_level; lev++) {

            int value = counts[t][lev];
            write_float(file, value);
        }
    }

    file.close();
}


std::vector<std::vector<int>> read_box_counts(std::string path,
                                              std::string in_file,
                                              int         max_time,
                                              int         min_time,
                                              int         min_level,
                                              int         max_level) {

    std::ifstream file = open_read(path, in_file);

    std::vector<std::vector<int>> read_counts;

    for (int t = 0; t < max_time - min_time + 1; t++) {
        std::vector<int> read_counts_t;

        for (int lev = min_level; lev <= max_level; lev++) {
            int value = static_cast<int>(read_float(file));
            read_counts_t.push_back(value);
        }

        read_counts.push_back(read_counts_t);
    }

    file.close();
    return read_counts;
}


double
calc_avg_rmse(std::vector<std::vector<std::vector<Box3D>>> const& actuals,
              std::vector<std::vector<std::vector<Box3D>>> const& regens) {

    std::vector<double> rmses;

    for (int t = 0; t < actuals.size(); t++) {

        auto const& actuals_t = actuals[t];
        auto const& regens_t  = regens[t];

        for (int l = 0; l < actuals_t.size(); l++) {

            auto const& actuals_l = actuals_t[l];
            auto const& regens_l  = regens_t[l];

            for (int box_idx = 0; box_idx < actuals_l.size(); box_idx++) {

                auto const& actual = actuals_l[box_idx];
                auto const& pred   = regens_l[box_idx];

                int xdim = actual.width();
                int ydim = actual.height();
                int zdim = actual.depth();

                double sum = 0.0;

                for (int k = 0; k < zdim; k++) {
                    for (int j = 0; j < ydim; j++) {
                        for (int i = 0; i < xdim; i++) {
                            double diff = actual(i, j, k) - pred(i, j, k);
                            sum += diff * diff;
                        }
                    }
                }

                double rmse = sqrt(sum / (xdim * ydim * zdim));
                rmses.push_back(rmse);
            }
        }
    }

    double mean_rmse =
        accumulate(rmses.begin(), rmses.end(), 0.0) / rmses.size();
    return mean_rmse;
}


void write_box_with_loc_dim(Volume3D           volume,
                            std::string        path,
                            std::string        out_file,
                            std::vector<float> location,
                            std::vector<float> dimension) {

    std::ofstream file = open_write(path, out_file);

    // Write location
    for (float coord : location) {
        write_float(file, coord);
    }

    // Write dimension
    for (float dim : dimension) {
        write_float(file, dim);
    }

    // Write volume data
    volume.iterate([&](float value, size_t, size_t, size_t) {
        write_float(file, value);
    });
    // TODO: Check and see if this can be used anywhere else

    file.close();
}


std::vector<std::pair<int, int16_t>> rle_encode(std::vector<bool>    mask,
                                                std::vector<int16_t> values) {

    std::vector<std::pair<int, int16_t>> rle;
    int                                  run_length = 0;
    int                                  val_idx    = 0;

    for (bool keep : mask) {

        if (keep) {

            rle.emplace_back(run_length, values[val_idx]);
            run_length = 0;
            val_idx++;

        } else {

            run_length++;
        }
    }

    return rle;
}


std::string serialize_compressed_wavelet(CompressedWavelet compressed) {
    std::string buffer;

    for (int dim : compressed.shape) {
        buffer.append(reinterpret_cast<const char*>(&dim), sizeof(int));
    }

    for (int dim : compressed.coeff_shape) {
        buffer.append(reinterpret_cast<const char*>(&dim), sizeof(int));
    }

    int rle_size = static_cast<int>(compressed.rle_encoded.size());
    buffer.append(reinterpret_cast<const char*>(&rle_size), sizeof(int));

    for (const auto& [run, val] : compressed.rle_encoded) {
        buffer.append(reinterpret_cast<const char*>(&run), sizeof(int));
        buffer.append(reinterpret_cast<const char*>(&val), sizeof(int16_t));
    }

    return buffer;
}


std::vector<double> wavelet_decompose(Box3D const& box) {
    int x = box.width();
    int y = box.height();
    int z = box.depth();

    Box3D               temp = box.clone();
    std::vector<double> flat;

    // Decompose along Z
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            std::vector<double> row(z);
            for (int k = 0; k < z; ++k)
                row[k] = temp(i, j, k);

            std::vector<double> low, high;
            for (size_t k = 0; k + 1 < row.size(); k += 2) {
                double a = row[k], b = row[k + 1];
                low.push_back((a + b) / 2.0);
                high.push_back((a - b) / 2.0);
            }

            size_t mid = low.size();
            for (size_t k = 0; k < mid; ++k)
                row[k] = low[k];
            for (size_t k = 0; k < high.size(); ++k)
                row[mid + k] = high[k];

            for (int k = 0; k < z; ++k)
                temp(i, j, k) = row[k];
        }
    }

    // Repeat for Y
    for (int i = 0; i < x; ++i) {
        for (int k = 0; k < z; ++k) {
            std::vector<double> col(y);
            for (int j = 0; j < y; ++j)
                col[j] = temp(i, j, k);

            std::vector<double> low, high;
            for (size_t j = 0; j + 1 < col.size(); j += 2) {
                double a = col[j], b = col[j + 1];
                low.push_back((a + b) / 2.0);
                high.push_back((a - b) / 2.0);
            }

            size_t mid = low.size();
            for (size_t j = 0; j < mid; ++j)
                col[j] = low[j];
            for (size_t j = 0; j < high.size(); ++j)
                col[mid + j] = high[j];

            for (int j = 0; j < y; ++j)
                temp(i, j, k) = col[j];
        }
    }

    // Repeat for X
    for (int j = 0; j < y; ++j) {
        for (int k = 0; k < z; ++k) {
            std::vector<double> depth(x);
            for (int i = 0; i < x; ++i)
                depth[i] = temp(i, j, k);

            std::vector<double> low, high;
            for (size_t i = 0; i + 1 < depth.size(); i += 2) {
                double a = depth[i], b = depth[i + 1];
                low.push_back((a + b) / 2.0);
                high.push_back((a - b) / 2.0);
            }

            size_t mid = low.size();
            for (size_t i = 0; i < mid; ++i)
                depth[i] = low[i];
            for (size_t i = 0; i < high.size(); ++i)
                depth[mid + i] = high[i];

            for (int i = 0; i < x; ++i)
                temp(i, j, k) = depth[i];
        }
    }

    // Flatten result
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            for (int k = 0; k < z; ++k)
                flat.push_back(temp(i, j, k));

    return flat;
}


CompressedWavelet wavelet_compress_rle(Box3D const& box,
                                       double       keep,
                                       // string wavelet,
                                       int         time,
                                       int         level,
                                       int         box_index,
                                       std::string compressed_dir) {

    // Wavelet decomposition
    std::vector<double> flat_coeffs = wavelet_decompose(box);

    // Thresholding
    double max_val =
        *std::max_element(flat_coeffs.begin(),
                          flat_coeffs.end(),
                          [](double a, double b) { return abs(a) < abs(b); });
    double thresh = max_val * (1 - keep);

    std::vector<bool>    mask;
    std::vector<int16_t> values;
    for (double val : flat_coeffs) {
        bool keep_val = abs(val) > thresh;
        mask.push_back(keep_val);
        if (keep_val) {
            values.push_back(static_cast<int16_t>(round(val * 1000.0) / 1000));
        }
    }

    // RLE encode
    std::vector<std::pair<int, int16_t>> rle_encoded = rle_encode(mask, values);

    // Create compressed structure
    CompressedWavelet compressed;
    compressed.shape       = { static_cast<int>(box.width()),
                               static_cast<int>(box.height()),
                               static_cast<int>(box.depth()) };
    compressed.coeff_shape = { static_cast<int>(flat_coeffs.size()) };
    compressed.rle_encoded = rle_encoded;

    // Serialize and compress using lzma
    std::string filename = compressed_dir + "compressed-wavelet-" +
                           std::to_string(time) + "-" + std::to_string(level) +
                           "-" + std::to_string(box_index) + ".xz";

    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        std::string serialized = serialize_compressed_wavelet(compressed);

        lzma_stream strm = LZMA_STREAM_INIT;
        lzma_ret    ret  = lzma_easy_encoder(
            &strm, 6 /* compression level */, LZMA_CHECK_CRC64);
        if (ret != LZMA_OK) {
            spdlog::error("Failed to initialize LZMA encoder: {}");
            exit(EXIT_FAILURE);
        }

        // Input
        strm.next_in  = reinterpret_cast<const uint8_t*>(serialized.data());
        strm.avail_in = serialized.size();

        // Output buffer (reserve slightly more than input)
        std::vector<uint8_t> outbuf(serialized.size() * 1.1 +
                                    128); // +128 safety
        strm.next_out  = outbuf.data();
        strm.avail_out = outbuf.size();

        // Compress
        ret = lzma_code(&strm, LZMA_FINISH);
        if (ret != LZMA_STREAM_END) {
            // This should be fatal
            spdlog::error("LZMA compression failed: {}");
            exit(EXIT_FAILURE);
        }

        lzma_end(&strm);

        // Write compressed data
        file.write(reinterpret_cast<const char*>(outbuf.data()),
                   outbuf.size() - strm.avail_out);
        file.close();
    }

    return compressed;
}


std::vector<float> rle_decode(std::vector<std::pair<int, int16_t>> rle_encoded,
                              int total_length) {

    std::vector<float> result(total_length, 0.0f);
    int                idx = 0;

    for (auto [run_length, value] : rle_encoded) {

        idx += run_length;
        if (idx < total_length) {
            result[idx] = value;
            ++idx;
        }
    }

    return result;
}


CompressedWavelet deserialize_compressed_wavelet(const std::string& data) {
    CompressedWavelet compressed;
    const char* ptr = data.data(); // Pointer to start of binary data

    auto read = [&](auto& out) {
        std::memcpy(&out, ptr, sizeof(out));
        ptr += sizeof(out);
    };

    // Shape (3 dimensions)
    for (int i = 0; i < 3; ++i) {
        int dim;
        std::memcpy(&dim, ptr, sizeof(dim));
        ptr += sizeof(dim);
        compressed.shape.push_back(dim);
    }

    // Coeff shape (1 dimension)
    for (int i = 0; i < 1; ++i) {
        int dim;
        std::memcpy(&dim, ptr, sizeof(dim));
        ptr += sizeof(dim);
        compressed.coeff_shape.push_back(dim);
    }

    // RLE size
    int rle_size;
    std::memcpy(&rle_size, ptr, sizeof(rle_size));
    ptr += sizeof(rle_size);
    compressed.rle_encoded.resize(rle_size);

    // RLE data
    for (int i = 0; i < rle_size; ++i) {
        int     run;
        int16_t val;
        std::memcpy(&run, ptr, sizeof(run));
        ptr += sizeof(run);
        std::memcpy(&val, ptr, sizeof(val));
        ptr += sizeof(val);
        compressed.rle_encoded[i] = { run, val };
    }

    return compressed;
}



Box3D inverse_wavelet_decompose(std::vector<float> flat, int x, int y, int z) {
    // Step 1: reshape the flat array into a 3D volume
    Box3D temp(x, y, z);

    int idx = 0;
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            for (int k = 0; k < z; ++k)
                temp(i, j, k) = flat[idx++];

    // Step 2: Inverse along X
    for (int j = 0; j < y; ++j) {
        for (int k = 0; k < z; ++k) {
            std::vector<double> depth(x);
            for (int i = 0; i < x; ++i)
                depth[i] = temp(i, j, k);

            std::vector<double> restored(x);
            int                 half = x / 2;
            for (int i = 0; i < half; ++i) {
                double avg          = depth[i];
                double diff         = depth[half + i];
                restored[2 * i]     = avg + diff;
                restored[2 * i + 1] = avg - diff;
            }

            for (int i = 0; i < x; ++i)
                temp(i, j, k) = restored[i];
        }
    }

    // Step 3: Inverse along Y
    for (int i = 0; i < x; ++i) {
        for (int k = 0; k < z; ++k) {
            std::vector<double> col(y);
            for (int j = 0; j < y; ++j)
                col[j] = temp(i, j, k);

            std::vector<double> restored(y);
            int                 half = y / 2;
            for (int j = 0; j < half; ++j) {
                double avg          = col[j];
                double diff         = col[half + j];
                restored[2 * j]     = avg + diff;
                restored[2 * j + 1] = avg - diff;
            }

            for (int j = 0; j < y; ++j)
                temp(i, j, k) = restored[j];
        }
    }

    // Step 4: Inverse along Z
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            std::vector<double> row(z);
            for (int k = 0; k < z; ++k)
                row[k] = temp(i, j, k);

            std::vector<double> restored(z);
            int                 half = z / 2;
            for (int k = 0; k < half; ++k) {
                double avg          = row[k];
                double diff         = row[half + k];
                restored[2 * k]     = avg + diff;
                restored[2 * k + 1] = avg - diff;
            }

            for (int k = 0; k < z; ++k)
                temp(i, j, k) = restored[k];
        }
    }

    return temp;
}


CompressedWavelet read_compressed_wavelet(const std::string& filename) {
    CompressedWavelet compressed;

    // Get file size and read contents
    std::error_code ec;
    auto            size = std::filesystem::file_size(filename, ec);
    if (ec) {
        spdlog::error("Error getting file size: {}", ec.message());
        exit(EXIT_FAILURE);
    }

    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        spdlog::error("Failed to open file: {}", filename);
        exit(EXIT_FAILURE);
    }

    std::vector<uint8_t> compressed_data(size);
    if (!file.read(reinterpret_cast<char*>(compressed_data.data()), size)) {
        spdlog::error("Failed to read file: ", filename);
        exit(EXIT_FAILURE);
    }

    // Initialize LZMA stream with RAII
    lzma_stream strm = LZMA_STREAM_INIT;
    if (lzma_stream_decoder(&strm, UINT64_MAX, LZMA_CONCATENATED) != LZMA_OK) {
        spdlog::error("Failed to initialize LZMA decoder.");
        exit(EXIT_FAILURE);
    }

    strm.next_in  = compressed_data.data();
    strm.avail_in = compressed_data.size();

    std::vector<uint8_t> decompressed_data;
    size_t               buffer_size = 4096;
    decompressed_data.resize(buffer_size);

    strm.next_out  = decompressed_data.data();
    strm.avail_out = buffer_size;

    while (true) {
        lzma_ret ret = lzma_code(&strm, LZMA_FINISH);
        if (ret == LZMA_STREAM_END) break;
        if (ret != LZMA_OK) {
            spdlog::error("LZMA decompression failed with code: {}", (int)ret);
            exit(EXIT_FAILURE);
        }

        // Buffer full, expand
        size_t old_size = decompressed_data.size();
        decompressed_data.resize(old_size * 2);
        strm.next_out  = decompressed_data.data() + old_size;
        strm.avail_out = old_size;
    }

    size_t decompressed_size = decompressed_data.size() - strm.avail_out;
    lzma_end(&strm);

    try {
        std::string serialized_data(
            reinterpret_cast<char*>(decompressed_data.data()),
            decompressed_size);
        compressed = deserialize_compressed_wavelet(serialized_data);
    } catch (const std::exception& e) {
        spdlog::error("Deserialization failed: {}", e.what());
        exit(EXIT_FAILURE);
    }

    return compressed;
}


int main(int argc, char* argv[]) {
    spdlog::set_level(spdlog::level::debug);

    if (argc != 7) {
        spdlog::info(
            "Usage: {} min_time max_time min_level max_level component keep",
            argv[0]);
        return EXIT_FAILURE;
    }

    // Parse command-line arguments
    int   min_time  = std::atoi(argv[1]);
    int   max_time  = std::atoi(argv[2]);
    int   min_level = std::atoi(argv[3]);
    int   max_level = std::atoi(argv[4]);
    int   component = std::atoi(argv[5]);
    float keep      = std::atof(argv[6]);

    std::string dir_prefix =
        "../../wholeVolumesNewFormat-" + std::to_string(component);

    std::filesystem::create_directory("../../wavelet/");


    spdlog::info("Processing data...");
    auto start          = std::chrono::high_resolution_clock::now();
    auto processed_data = process_data(
        dir_prefix, min_time, max_time, min_level, max_level, component);

    auto& boxes      = processed_data.boxes;
    auto& locations  = processed_data.locations;
    auto& dimensions = processed_data.dimensions;
    auto& box_counts = processed_data.box_counts;
    auto& min_value  = processed_data.min_value;
    auto& max_value  = processed_data.max_value;

    write_loc_dim_to_bin(locations,
                         compressed_dir,
                         "locations.raw",
                         min_time,
                         max_time,
                         min_level,
                         max_level);
    write_loc_dim_to_bin(dimensions,
                         compressed_dir,
                         "dimensions.raw",
                         min_time,
                         max_time,
                         min_level,
                         max_level);
    write_box_counts(box_counts,
                     compressed_dir,
                     "boxcounts.raw",
                     min_time,
                     max_time,
                     min_level,
                     max_level);

    auto end      = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end - start).count();
    spdlog::info(
        "Successfully processed data in {} seconds. Beginning compression...",
        duration);

    auto start1 = std::chrono::high_resolution_clock::now();
    for (int t = 0; t <= max_time - min_time; t++) {
        for (int lev = 0; lev <= max_level - min_level; lev++) {
            for (int box_idx = 0; box_idx < boxes[t][lev].size(); box_idx++) {

                Box3D const&      current_box = boxes[t][lev][box_idx];
                CompressedWavelet compressed  = wavelet_compress_rle(
                    current_box, keep, t, lev, box_idx, compressed_dir);
            }
        }
    }

    auto end1      = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration<double>(end1 - start1).count();
    spdlog::info(
        "Compression completed in {} seconds. Beginning decompression...",
        duration1);

    auto start2      = std::chrono::high_resolution_clock::now();
    auto counts_read = read_box_counts(compressed_dir,
                                       "boxcounts.raw",
                                       min_time,
                                       max_time,
                                       min_level,
                                       max_level);

    std::vector<std::vector<std::vector<Box3D>>> regen_boxes;

    for (int t = 0; t <= max_time - min_time; t++) {
        std::vector<std::vector<Box3D>> regen_boxes_t;
        for (int lev = 0; lev <= max_level - min_level; lev++) {
            std::vector<Box3D> regen_boxes_l;
            for (int box_idx = 0; box_idx < boxes[t][lev].size(); box_idx++) {

                std::string file_path = compressed_dir + "compressed-wavelet-" +
                                        std::to_string(t) + "-" +
                                        std::to_string(lev) + "-" +
                                        std::to_string(box_idx) + ".xz";

                CompressedWavelet read_compressed =
                    read_compressed_wavelet(file_path);

                auto read_flat = rle_decode(read_compressed.rle_encoded,
                                            read_compressed.coeff_shape[0]);

                auto const& dims      = read_compressed.shape;
                Box3D       regen_box = inverse_wavelet_decompose(
                    read_flat, dims[0], dims[1], dims[2]);

                regen_boxes_l.push_back(std::move(regen_box));
            }
            regen_boxes_t.push_back(std::move(regen_boxes_l));
        }
        regen_boxes.push_back(std::move(regen_boxes_t));
    }

    auto end2      = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration<double>(end2 - start2).count();
    spdlog::info("Decompression completed in {} seconds. Calculating loss...",
                 duration2);

    double rmse = calc_avg_rmse(boxes, regen_boxes);
    spdlog::info("RMSE = {}", rmse);

    double range = max_value - min_value;
    double loss = rmse / range;
    spdlog::info("Adjusted loss = {}", loss);
}


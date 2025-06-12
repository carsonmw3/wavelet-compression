#include "compressor.h"

#include <lzma.h>
#include <spdlog/spdlog.h>

#include <fstream>
#include <string>

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


CompressedWavelet compress(Box3D const& box,
                           double       keep,
                           int          time,
                           int          level,
                           int          box_index,
                           std::string  compressed_dir) {

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

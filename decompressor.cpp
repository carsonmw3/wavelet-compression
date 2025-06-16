#include "decompressor.h"

#include <lzma.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <fstream>
#include <string>

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


static CompressedWavelet deserialize_compressed_wavelet(const std::string& data) {
    CompressedWavelet compressed;
    const char*       ptr = data.data(); // Pointer to start of binary data

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

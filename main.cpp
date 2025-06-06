#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <tuple>
#include <numeric>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <cstring>
#include <lzma.h>
#include <filesystem>
#include <chrono>
#include <algorithm>
#include <cstdlib>


using namespace std;


// Hyperparameters
// int min_time = 74;
// int max_time = 74;
// int min_level = 0;
// int max_level = 3;
// int component = 6;

string compressed_dir = "../../wavelet/";

// string wavelet = "db3";
// float keep = 0.9999;



// Typedef for readability
typedef vector<vector<vector<float>>> Box3D;
typedef vector<vector<vector<float>>> Volume3D;



struct CompressedWavelet {
    vector<int> shape;
    vector<int> coeff_shape;
    // vector<vector<int>> coeff_slices;
    vector<pair<int, int16_t>> rle_encoded;
    // string wavelet;
};



// Function to create an empty 3D block filled with a junk value
Box3D create_block(int block_dim, float junk) {
    return Box3D(block_dim, vector<vector<float>>(block_dim, vector<float>(block_dim, junk)));
}


// Function to split the box into blocks
vector<Box3D> to_blocks(const Box3D& box, int block_dim, float junk) {

    int xdim = box.size();
    int ydim = box[0].size();
    int zdim = box[0][0].size();

    int x_steps = (xdim + block_dim - 1) / block_dim;
    int y_steps = (ydim + block_dim - 1) / block_dim;
    int z_steps = (zdim + block_dim - 1) / block_dim;

    vector<Box3D> blocks;

    for (int k = 0; k < z_steps; ++k) {
        for (int j = 0; j < y_steps; ++j) {
            for (int i = 0; i < x_steps; ++i) {

                int x_lo = i * block_dim;
                int y_lo = j * block_dim;
                int z_lo = k * block_dim;

                int x_hi = min(x_lo + block_dim, xdim);
                int y_hi = min(y_lo + block_dim, ydim);
                int z_hi = min(z_lo + block_dim, zdim);

                Box3D block = create_block(block_dim, junk);

                for (int x = x_lo; x < x_hi; ++x) {
                    for (int y = y_lo; y < y_hi; ++y) {
                        for (int z = z_lo; z < z_hi; ++z) {
                            block[x - x_lo][y - y_lo][z - z_lo] = box[x][y][z];
                        }
                    }
                }

                blocks.push_back(block);

            }
        }
    }

    return blocks;
}



// Read a float from file
float read_float(std::ifstream& file) {
    float val;
    file.read(reinterpret_cast<char*>(&val), sizeof(float));
    return val;
}



tuple<
    vector<vector<vector<Box3D>>>,               // boxes[time][level][box]
    vector<vector<vector<vector<int>>>>, // locations[time][level][box]
    vector<vector<vector<vector<int>>>>, // dimensions[time][level][box]
    vector<vector<int>>                          // box_counts[time][level]
    >
process_data(const string& dir_prefix, int min_time, int max_time, int min_level, int max_level, int component) {

    vector<vector<vector<Box3D>>> boxes;
    vector<vector<vector<vector<int>>>> locations;
    vector<vector<vector<vector<int>>>> dimensions;
    vector<vector<int>> box_counts;

    for (int t = min_time; t <= max_time; ++t) {
        vector<vector<Box3D>> boxes_t;
        vector<vector<vector<int>>> locations_t;
        vector<vector<vector<int>>> dimensions_t;
        vector<int> box_counts_t;

        for (int l = min_level; l <= max_level; ++l) {
            string filename = dir_prefix + "-" + to_string(l) + "/" +
                              to_string(t) + "-wholeNewFormat-" +
                              to_string(component) + "-" + to_string(l) + ".raw";

            ifstream file(filename, ios::binary);
            if (!file.is_open()) {
                cerr << "Could not open file: " << filename << endl;
                continue;
            }

            int box_count = 0;
            vector<Box3D> boxes_l;
            vector<vector<int>> locations_l;
            vector<vector<int>> dimensions_l;

            while (file.peek() != EOF) {
                float x_f;
                file.read(reinterpret_cast<char*>(&x_f), sizeof(float));
                if (file.gcount() < sizeof(float)) break;
                float y_f = read_float(file);
                float z_f = read_float(file);

                vector<int> current_loc;
                int x = static_cast<int>(x_f);
                current_loc.push_back(x);
                int y = static_cast<int>(y_f);
                current_loc.push_back(y);
                int z = static_cast<int>(z_f);
                current_loc.push_back(z);
                locations_l.push_back(current_loc);

                vector<int> current_dim;
                int xdim = static_cast<int>(read_float(file));
                current_dim.push_back(x);
                int ydim = static_cast<int>(read_float(file));
                current_dim.push_back(y);
                int zdim = static_cast<int>(read_float(file));
                current_dim.push_back(z);
                dimensions_l.push_back(current_dim);

                Box3D box(xdim, vector<vector<float>>(ydim, vector<float>(zdim)));

                for (int k = 0; k < zdim; ++k) {
                    for (int j = 0; j < ydim; ++j) {
                        for (int i = 0; i < xdim; ++i) {
                            float value = read_float(file);
                            box[i][j][k] = value;
                        }
                    }
                }

                boxes_l.push_back(box);
                ++box_count;
            }

            boxes_t.push_back(boxes_l);
            locations_t.push_back(locations_l);
            dimensions_t.push_back(dimensions_l);
            box_counts_t.push_back(box_count);

        }

        boxes.push_back(boxes_t);
        locations.push_back(locations_t);
        dimensions.push_back(dimensions_t);
        box_counts.push_back(box_counts_t);

    }

    return make_tuple(boxes, locations, dimensions, box_counts);

}





void write_loc_dim_to_bin(vector<vector<vector<vector<int>>>> data, string path, string out_file,
                          int max_time, int min_time, int min_level, int max_level) {

    string filename = path + out_file;

    ofstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Could not open file: " << filename << endl;
    }

    for (int t = 0; t < max_time-min_time+1; t++) {

        for (int lev = min_level; lev <= max_level; lev++) {

            for (int box_idx = 0; box_idx < data[t][lev].size(); box_idx++) {

                for (int coord = 0; coord < 3; coord++) {

                    float value = data[t][lev][box_idx][coord];
                    file.write(reinterpret_cast<char*>(&value), sizeof(float));

                }
            }
        }
    }

    file.close();
}



vector<vector<vector<vector<int>>>> read_loc_dim_from_bin(string path, string in_file, vector<vector<int>> counts,
                                                          int max_time, int min_time, int min_level, int max_level) {

    string filename = path + in_file;
    ifstream file(filename, ios::binary);

    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return {};
    }

    vector<vector<vector<vector<int>>>> out;

    for (int t = 0; t < max_time-min_time+1; t++) {
        vector<vector<vector<int>>> out_t;

        for (int lev = min_level; lev <= max_level; lev++) {
            vector<vector<int>> out_l;
            int curr_num_boxes = counts[t][lev];

            for (int box_idx = 0; box_idx < curr_num_boxes; box_idx++) {
                vector<int> current_coords;

                for (int coord = 0; coord < 3; coord++) {
                    float value;
                    file.read(reinterpret_cast<char*>(&value), sizeof(float));
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




void write_box_counts(vector<vector<int>> counts, string path, string out_file,
                      int max_time, int min_time, int min_level, int max_level) {

    string filename = path + out_file;
    ofstream file(filename, ios::binary);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    for (int t = 0; t < max_time-min_time+1; t++) {

        for (int lev = min_level; lev <= max_level; lev++) {

            int value = counts[t][lev];
            file.write(reinterpret_cast<const char*>(&value), sizeof(float));

        }
    }

    file.close();
}



vector<vector<int>> read_box_counts(string path, string in_file,
                                    int max_time, int min_time, int min_level, int max_level) {

    string filename = path + in_file;
    ifstream file(filename, ios::binary);

    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return {};
    }

    vector<vector<int>> read_counts;

    for (int t = 0; t < max_time-min_time+1; t++) {
        vector<int> read_counts_t;

        for (int lev = min_level; lev <= max_level; lev++) {
            int value;
            file.read(reinterpret_cast<char*>(&value), sizeof(float));
            read_counts_t.push_back(value);
        }

        read_counts.push_back(read_counts_t);
    }

    file.close();
    return read_counts;
}





double calc_avg_rmse(vector<vector<vector<Box3D>>> actuals, vector<vector<vector<Box3D>>> regens) {

    vector<double> rmses;

    for (int t = 0; t < actuals.size(); t++) {

        auto actuals_t = actuals[t];
        auto regens_t = regens[t];

        for (int l = 0; l < actuals_t.size(); l++) {

            auto actuals_l = actuals_t[l];
            auto regens_l = regens_t[l];

            for (int box_idx = 0; box_idx < actuals_l.size(); box_idx++) {

                auto actual = actuals_l[box_idx];
                auto pred = regens_l[box_idx];

                int xdim = actual.size();
                int ydim = actual[0].size();
                int zdim = actual[0][0].size();

                double sum = 0.0;

                for (int k = 0; k < zdim; k++) {
                    for (int j = 0; j < ydim; j++) {
                        for (int i = 0; i < xdim; i++) {
                            double diff = actual[i][j][k] - pred[i][j][k];
                            sum += diff * diff;
                        }
                    }
                }

                double rmse = sqrt(sum / (xdim * ydim * zdim));
                rmses.push_back(rmse);

            }

        }

    }

    double mean_rmse = accumulate(rmses.begin(), rmses.end(), 0.0) / rmses.size();
    return mean_rmse;

}




void write_box_with_loc_dim(Volume3D volume, string path, string out_file, vector<float> location, vector<float> dimension) {

    string filename = path + out_file;
    ofstream file(filename, ios::binary);

    if (!file) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write location
    for (float coord : location) {
        file.write(reinterpret_cast<const char*>(&coord), sizeof(float));
    }

    // Write dimension
    for (float dim : dimension) {
        file.write(reinterpret_cast<const char*>(&dim), sizeof(float));
    }

    // Write volume data
    int xdim = volume.size();
    int ydim = volume[0].size();
    int zdim = volume[0][0].size();

    for (int z = 0; z < zdim; z++) {
        for (int y = 0; y < ydim; y++) {
            for (int x = 0; x < xdim; x++) {
                float val = volume[x][y][z];
                file.write(reinterpret_cast<const char*>(&val), sizeof(float));
            }
        }
    }

    file.close();
}




vector<pair<int, int16_t>> rle_encode(vector<bool> mask, vector<int16_t> values) {

    vector<pair<int, int16_t>> rle;
    int run_length = 0;
    int val_idx = 0;

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



string serialize_compressed_wavelet(CompressedWavelet compressed) {

    ostringstream oss(ios::binary);

    // Shape
    for (int dim : compressed.shape) {
        oss.write(reinterpret_cast<const char*>(&dim), sizeof(int));
    }

    // Coeff shape
    for (int dim : compressed.coeff_shape) {
        oss.write(reinterpret_cast<const char*>(&dim), sizeof(int));
    }

    // RLE data
    int rle_size = static_cast<int>(compressed.rle_encoded.size());
    oss.write(reinterpret_cast<const char*>(&rle_size), sizeof(int));
    for (const auto& [run, val] : compressed.rle_encoded) {
        oss.write(reinterpret_cast<const char*>(&run), sizeof(int));
        oss.write(reinterpret_cast<const char*>(&val), sizeof(int16_t));
    }

    // Wavelet name
    // int str_len = static_cast<int>(compressed.wavelet.size());
    // oss.write(reinterpret_cast<const char*>(&str_len), sizeof(int));
    // oss.write(compressed.wavelet.data(), str_len);

    return oss.str();

}



vector<double> wavelet_decompose(const Box3D& box) {
    int x = box.size();
    int y = box[0].size();
    int z = box[0][0].size();

    Box3D temp = box;
    vector<double> flat;

    // Decompose along Z
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            vector<double> row(z);
            for (int k = 0; k < z; ++k)
                row[k] = temp[i][j][k];

            vector<double> low, high;
            for (size_t k = 0; k + 1 < row.size(); k += 2) {
                double a = row[k], b = row[k+1];
                low.push_back((a + b) / 2.0);
                high.push_back((a - b) / 2.0);
            }

            size_t mid = low.size();
            for (size_t k = 0; k < mid; ++k)
                row[k] = low[k];
            for (size_t k = 0; k < high.size(); ++k)
                row[mid + k] = high[k];

            for (int k = 0; k < z; ++k)
                temp[i][j][k] = row[k];
        }
    }

    // Repeat for Y
    for (int i = 0; i < x; ++i) {
        for (int k = 0; k < z; ++k) {
            vector<double> col(y);
            for (int j = 0; j < y; ++j)
                col[j] = temp[i][j][k];

            vector<double> low, high;
            for (size_t j = 0; j + 1 < col.size(); j += 2) {
                double a = col[j], b = col[j+1];
                low.push_back((a + b) / 2.0);
                high.push_back((a - b) / 2.0);
            }

            size_t mid = low.size();
            for (size_t j = 0; j < mid; ++j)
                col[j] = low[j];
            for (size_t j = 0; j < high.size(); ++j)
                col[mid + j] = high[j];

            for (int j = 0; j < y; ++j)
                temp[i][j][k] = col[j];
        }
    }

    // Repeat for X
    for (int j = 0; j < y; ++j) {
        for (int k = 0; k < z; ++k) {
            vector<double> depth(x);
            for (int i = 0; i < x; ++i)
                depth[i] = temp[i][j][k];

            vector<double> low, high;
            for (size_t i = 0; i + 1 < depth.size(); i += 2) {
                double a = depth[i], b = depth[i+1];
                low.push_back((a + b) / 2.0);
                high.push_back((a - b) / 2.0);
            }

            size_t mid = low.size();
            for (size_t i = 0; i < mid; ++i)
                depth[i] = low[i];
            for (size_t i = 0; i < high.size(); ++i)
                depth[mid + i] = high[i];

            for (int i = 0; i < x; ++i)
                temp[i][j][k] = depth[i];
        }
    }

    // Flatten result
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            for (int k = 0; k < z; ++k)
                flat.push_back(temp[i][j][k]);

    return flat;
}




CompressedWavelet wavelet_compress_rle(
    Box3D box,
    double keep,
    // string wavelet,
    int time,
    int level,
    int box_index,
    string compressed_dir
    ) {

    // Wavelet decomposition
    vector<double> flat_coeffs = wavelet_decompose(box);

    // Thresholding
    double max_val = *std::max_element(flat_coeffs.begin(), flat_coeffs.end(),
                                       [](double a, double b) { return abs(a) < abs(b); });
    double thresh = max_val * (1 - keep);

    vector<bool> mask;
    vector<int16_t> values;
    for (double val : flat_coeffs) {
        bool keep_val = abs(val) > thresh;
        mask.push_back(keep_val);
        if (keep_val) {
            values.push_back(static_cast<int16_t>(round(val*1000.0) / 1000));
        }
    }

    // RLE encode
    vector<pair<int, int16_t>> rle_encoded = rle_encode(mask, values);

    // Create compressed structure
    CompressedWavelet compressed;
    compressed.shape = { static_cast<int>(box.size()), static_cast<int>(box[0].size()), static_cast<int>(box[0][0].size()) };
    compressed.coeff_shape = { static_cast<int>(flat_coeffs.size()) };
    // compressed.coeff_slices = {}; // Not implemented
    compressed.rle_encoded = rle_encoded;
    // compressed.wavelet = wavelet;

    // Serialize and compress using lzma
    string filename = compressed_dir + "compressed-wavelet-" + to_string(time) + "-" +
                      to_string(level) + "-" + to_string(box_index) + ".xz";

    ofstream file(filename, ios::binary);
    if (file.is_open()) {
        string serialized = serialize_compressed_wavelet(compressed);

        lzma_stream strm = LZMA_STREAM_INIT;
        lzma_ret ret = lzma_easy_encoder(&strm, 6 /* compression level */, LZMA_CHECK_CRC64);
        if (ret != LZMA_OK) {
            cerr << "Failed to initialize LZMA encoder: " << ret << endl;
            return compressed;
        }

        // Input
        strm.next_in = reinterpret_cast<const uint8_t*>(serialized.data());
        strm.avail_in = serialized.size();

        // Output buffer (reserve slightly more than input)
        vector<uint8_t> outbuf(serialized.size() * 1.1 + 128);  // +128 safety
        strm.next_out = outbuf.data();
        strm.avail_out = outbuf.size();

        // Compress
        ret = lzma_code(&strm, LZMA_FINISH);
        if (ret != LZMA_STREAM_END) {
            cerr << "LZMA compression failed: " << ret << endl;
            lzma_end(&strm);
            return compressed;
        }

        lzma_end(&strm);

        // Write compressed data
        file.write(reinterpret_cast<const char*>(outbuf.data()), outbuf.size() - strm.avail_out);
        file.close();
    }

    return compressed;
}




vector<float> rle_decode(vector<pair<int, int16_t>> rle_encoded, int total_length) {

    vector<float> result(total_length, 0.0f);
    int idx = 0;

    for (auto [run_length, value] : rle_encoded) {

        idx += run_length;
        if (idx < total_length) {
            result[idx] = value;
            ++idx;
        }
    }

    return result;
}



CompressedWavelet deserialize_compressed_wavelet(const string& data) {
    CompressedWavelet compressed;
    istringstream iss(data, ios::binary);

    // Shape (expecting 3 dimensions)
    for (int i = 0; i < 3; ++i) {
        int dim;
        iss.read(reinterpret_cast<char*>(&dim), sizeof(int));
        compressed.shape.push_back(dim);
    }

    // Coeff shape (expecting 1 dimension)
    for (int i = 0; i < 1; ++i) {
        int dim;
        iss.read(reinterpret_cast<char*>(&dim), sizeof(int));
        compressed.coeff_shape.push_back(dim);
    }

    // RLE data
    int rle_size;
    iss.read(reinterpret_cast<char*>(&rle_size), sizeof(int));
    compressed.rle_encoded.resize(rle_size);

    for (int i = 0; i < rle_size; ++i) {
        int run;
        int16_t val;
        iss.read(reinterpret_cast<char*>(&run), sizeof(int));
        iss.read(reinterpret_cast<char*>(&val), sizeof(int16_t));
        compressed.rle_encoded[i] = {run, val};
    }

    // Wavelet name
    // int str_len;
    // iss.read(reinterpret_cast<char*>(&str_len), sizeof(int));
    // string wavelet(str_len, '\0');
    // iss.read(&wavelet[0], str_len);
    // compressed.wavelet = wavelet;

    return compressed;
}



Box3D inverse_wavelet_decompose(vector<float> flat, int x, int y, int z) {
    // Step 1: reshape the flat array into a 3D volume
    vector<vector<vector<float>>> box(x, vector<vector<float>>(y, vector<float>(z)));
    Box3D temp = box;
    int idx = 0;
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            for (int k = 0; k < z; ++k)
                temp[i][j][k] = flat[idx++];

    // Step 2: Inverse along X
    for (int j = 0; j < y; ++j) {
        for (int k = 0; k < z; ++k) {
            vector<double> depth(x);
            for (int i = 0; i < x; ++i)
                depth[i] = temp[i][j][k];

            vector<double> restored(x);
            int half = x / 2;
            for (int i = 0; i < half; ++i) {
                double avg = depth[i];
                double diff = depth[half + i];
                restored[2*i]     = avg + diff;
                restored[2*i + 1] = avg - diff;
            }

            for (int i = 0; i < x; ++i)
                temp[i][j][k] = restored[i];
        }
    }

    // Step 3: Inverse along Y
    for (int i = 0; i < x; ++i) {
        for (int k = 0; k < z; ++k) {
            vector<double> col(y);
            for (int j = 0; j < y; ++j)
                col[j] = temp[i][j][k];

            vector<double> restored(y);
            int half = y / 2;
            for (int j = 0; j < half; ++j) {
                double avg = col[j];
                double diff = col[half + j];
                restored[2*j]     = avg + diff;
                restored[2*j + 1] = avg - diff;
            }

            for (int j = 0; j < y; ++j)
                temp[i][j][k] = restored[j];
        }
    }

    // Step 4: Inverse along Z
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            vector<double> row(z);
            for (int k = 0; k < z; ++k)
                row[k] = temp[i][j][k];

            vector<double> restored(z);
            int half = z / 2;
            for (int k = 0; k < half; ++k) {
                double avg = row[k];
                double diff = row[half + k];
                restored[2*k]     = avg + diff;
                restored[2*k + 1] = avg - diff;
            }

            for (int k = 0; k < z; ++k)
                temp[i][j][k] = restored[k];
        }
    }

    return temp;
}




// CompressedWavelet read_compressed_wavelet(const string& filename) {
//     CompressedWavelet compressed;

//     // Open file and read into buffer
//     ifstream file(filename, ios::binary | ios::ate);
//     if (!file.is_open()) {
//         cerr << "Failed to open file: " << filename << endl;
//         return compressed;
//     }

//     streamsize size = file.tellg();
//     file.seekg(0, ios::beg);

//     vector<uint8_t> compressed_data(size);
//     if (!file.read(reinterpret_cast<char*>(compressed_data.data()), size)) {
//         cerr << "Failed to read compressed data from file: " << filename << endl;
//         return compressed;
//     }

//     file.close();

//     // Initialize LZMA stream
//     lzma_stream strm = LZMA_STREAM_INIT;
//     if (lzma_stream_decoder(&strm, UINT64_MAX, LZMA_CONCATENATED) != LZMA_OK) {
//         cerr << "Failed to initialize LZMA decoder." << endl;
//         return compressed;
//     }

//     strm.next_in = compressed_data.data();
//     strm.avail_in = compressed_data.size();

//     // Set a max expected decompressed size to prevent infinite buffer growth (adjust as needed)
//     const size_t max_decompressed_size = 50 * 1024 * 1024; // 50 MB
//     vector<uint8_t> decompressed_data(max_decompressed_size);

//     strm.next_out = decompressed_data.data();
//     strm.avail_out = decompressed_data.size();

//     lzma_ret ret = lzma_code(&strm, LZMA_FINISH);
//     if (ret != LZMA_STREAM_END) {
//         cerr << "LZMA decompression failed (file: " << filename << "), code: " << ret << endl;
//         lzma_end(&strm);
//         return compressed;
//     }

//     size_t decompressed_size = decompressed_data.size() - strm.avail_out;
//     lzma_end(&strm);

//     // Deserialize
//     string serialized_data(reinterpret_cast<char*>(decompressed_data.data()), decompressed_size);
//     try {
//         compressed = deserialize_compressed_wavelet(serialized_data);
//     } catch (const std::exception& e) {
//         cerr << "Deserialization failed for file: " << filename << " with error: " << e.what() << endl;
//     }

//     return compressed;
// }



CompressedWavelet read_compressed_wavelet(const std::string& filename) {
    CompressedWavelet compressed;

    // Get file size and read contents
    std::error_code ec;
    auto size = std::filesystem::file_size(filename, ec);
    if (ec) {
        std::cerr << "Error getting file size: " << ec.message() << std::endl;
        return compressed;
    }

    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return compressed;
    }

    std::vector<uint8_t> compressed_data(size);
    if (!file.read(reinterpret_cast<char*>(compressed_data.data()), size)) {
        std::cerr << "Failed to read file: " << filename << std::endl;
        return compressed;
    }

    // Initialize LZMA stream with RAII
    lzma_stream strm = LZMA_STREAM_INIT;
    if (lzma_stream_decoder(&strm, UINT64_MAX, LZMA_CONCATENATED) != LZMA_OK) {
        std::cerr << "Failed to initialize LZMA decoder." << std::endl;
        return compressed;
    }

    strm.next_in = compressed_data.data();
    strm.avail_in = compressed_data.size();

    std::vector<uint8_t> decompressed_data;
    size_t buffer_size = 4096;
    decompressed_data.resize(buffer_size);

    strm.next_out = decompressed_data.data();
    strm.avail_out = buffer_size;

    while (true) {
        lzma_ret ret = lzma_code(&strm, LZMA_FINISH);
        if (ret == LZMA_STREAM_END) break;
        if (ret != LZMA_OK) {
            std::cerr << "LZMA decompression failed with code: " << ret << std::endl;
            lzma_end(&strm);
            return compressed;
        }

        // Buffer full, expand
        size_t old_size = decompressed_data.size();
        decompressed_data.resize(old_size * 2);
        strm.next_out = decompressed_data.data() + old_size;
        strm.avail_out = old_size;
    }

    size_t decompressed_size = decompressed_data.size() - strm.avail_out;
    lzma_end(&strm);

    try {
        std::string serialized_data(reinterpret_cast<char*>(decompressed_data.data()), decompressed_size);
        compressed = deserialize_compressed_wavelet(serialized_data);
    } catch (const std::exception& e) {
        std::cerr << "Deserialization failed: " << e.what() << std::endl;
    }

    return compressed;
}







int main (int argc, char* argv[]) {

    if (argc != 7) {
        cerr << "Usage: " << argv[0] << " min_time max_time min_level max_level component keep\n";
        return 1;
    }

    // Parse command-line arguments
    int min_time = std::atoi(argv[1]);
    int max_time = std::atoi(argv[2]);
    int min_level = std::atoi(argv[3]);
    int max_level = std::atoi(argv[4]);
    int component = std::atoi(argv[5]);
    float keep = std::atof(argv[6]);

    string dir_prefix = "../../wholeVolumesNewFormat-" + to_string(component);

    filesystem::create_directory("../../wavelet/");

    tuple<
        vector<vector<vector<Box3D>>>,
        vector<vector<vector<vector<int>>>>,
        vector<vector<vector<vector<int>>>>,
        vector<vector<int>>
        > processed_data;

    cout << "Processing data..." << endl;
    auto start = chrono::high_resolution_clock::now();
    processed_data = process_data(dir_prefix, min_time, max_time, min_level, max_level, component);

    vector<vector<vector<Box3D>>> boxes = get<0>(processed_data);
    vector<vector<vector<vector<int>>>> locations = get<1>(processed_data);
    vector<vector<vector<vector<int>>>> dimensions = get<2>(processed_data);
    vector<vector<int>> box_counts = get<3>(processed_data);

    write_loc_dim_to_bin(locations, compressed_dir, "locations.raw", min_time, max_time, min_level, max_level);
    write_loc_dim_to_bin(dimensions, compressed_dir, "dimensions.raw", min_time, max_time, min_level, max_level);
    write_box_counts(box_counts, compressed_dir, "boxcounts.raw", min_time, max_time, min_level, max_level);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double>(end - start).count();
    cout << "Successfully processed data in " << to_string(duration) << " seconds. Beginning compression..." << endl;

    auto start1 = chrono::high_resolution_clock::now();
    for (int t = 0; t <= max_time-min_time; t++) {
        for (int lev = 0; lev <= max_level-min_level; lev++) {
            for (int box_idx = 0; box_idx < boxes[t][lev].size(); box_idx++) {

                Box3D current_box = boxes[t][lev][box_idx];
                CompressedWavelet compressed = wavelet_compress_rle(current_box, keep, t, lev, box_idx, compressed_dir);

            }
        }
    }

    auto end1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration<double>(end1 - start1).count();
    cout << "Compression completed in " << to_string(duration1) << " seconds. Beginning decompression..." << endl;

    auto start2 = chrono::high_resolution_clock::now();
    vector<vector<int>> counts_read = read_box_counts(compressed_dir, "boxcounts.raw", min_time, max_time, min_level, max_level);
    // vector<vector<vector<vector<int>>>> dims_read = read_loc_dim_from_bin(compressed_dir, "dimensions.raw", counts_read);
    vector<vector<vector<Box3D>>> regen_boxes;

    for (int t = 0; t <= max_time-min_time; t++) {
        vector<vector<Box3D>> regen_boxes_t;
        for (int lev = 0; lev <= max_level-min_level; lev++) {
            vector<Box3D> regen_boxes_l;
            for (int box_idx = 0; box_idx < boxes[t][lev].size(); box_idx++) {

                string file_path = compressed_dir + "compressed-wavelet-" + to_string(t) + "-" +
                                   to_string(lev) + "-" + to_string(box_idx) + ".xz";

                // auto startRead = chrono::high_resolution_clock::now();
                CompressedWavelet read_compressed = read_compressed_wavelet(file_path);
                // auto endRead = chrono::high_resolution_clock::now();
                // auto durationRead = chrono::duration<double>(endRead - startRead).count();
                // cout << "Read: " << to_string(durationRead) << endl;

                // auto startDecode = chrono::high_resolution_clock::now();
                vector<float> read_flat = rle_decode(read_compressed.rle_encoded, read_compressed.coeff_shape[0]);
                // auto endDecode = chrono::high_resolution_clock::now();
                // auto durationDecode = chrono::duration<double>(endDecode - startDecode).count();
                // cout << "Decode: " << to_string(durationDecode) << endl;

                // auto startWavelet = chrono::high_resolution_clock::now();
                vector<int> dims = read_compressed.shape;
                Box3D regen_box = inverse_wavelet_decompose(read_flat, dims[0], dims[1], dims[2]);
                // auto endWavelet = chrono::high_resolution_clock::now();
                // auto durationWavelet = chrono::duration<double>(endWavelet - startWavelet).count();
                // cout << "Wavelet: " << to_string(durationWavelet) << endl;

                regen_boxes_l.push_back(regen_box);
            }
            regen_boxes_t.push_back(regen_boxes_l);
        }
        regen_boxes.push_back(regen_boxes_t);
    }

    auto end2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration<double>(end2 - start2).count();
    cout << "Decompression completed in " << to_string(duration2) << " seconds. Calculating loss..." << endl;

    double rmse = calc_avg_rmse(boxes, regen_boxes);
    cout << "RMSE = " << rmse << endl;

}




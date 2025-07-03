#pragma once

#include <filesystem>
#include <string>
#include <unistd.h>

class TempDir {
    std::filesystem::path m_path;

    static std::string generate_temp_dir() {
        std::string name_template =
            std::filesystem::temp_directory_path() / "wavelet_compress.XXXXXX";
        return mkdtemp(name_template.data());
    }

public:
    TempDir() : m_path(generate_temp_dir()) { }

    // Permit no moves nor copy
    TempDir(const TempDir&)            = delete;
    TempDir& operator=(const TempDir&) = delete;
    TempDir(TempDir&&)                 = delete;
    TempDir& operator=(TempDir&&)      = delete;

    std::filesystem::path const& path() const { return m_path; }

    ~TempDir() {
        if (!m_path.empty() and std::filesystem::is_directory(m_path)) {
            try {
                std::filesystem::remove_all(m_path);
            } catch (...) {
                // Do not throw in destructor
            }
        }
    }
};

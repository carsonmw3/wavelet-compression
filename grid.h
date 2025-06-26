#pragma once

#include <vector>
#include <stdexcept>

template <class T>
class Grid3D {
private:
    size_t         m_width  = 0;
    size_t         m_height = 0;
    size_t         m_depth  = 0;
    std::vector<T> m_data;

    // Convert 3D indices to a flat 1D index
    size_t index(size_t x, size_t y, size_t z) const {
        if (x >= m_width || y >= m_height || z >= m_depth)
            throw std::out_of_range("Grid3D index out of range");
        return x + m_width * (y + m_height * z);
    }

public:
    // Constructors
    Grid3D() { }
    Grid3D(size_t width, size_t height, size_t depth, T initial_value = T())
        : m_width(width),
          m_height(height),
          m_depth(depth),
          m_data(width * height * depth, initial_value) { }

    // Disable expensive copying
    Grid3D(Grid3D const&)            = delete;
    Grid3D& operator=(Grid3D const&) = delete;

    // Enable moving
    Grid3D(Grid3D&&)            = default;
    Grid3D& operator=(Grid3D&&) = default;

    // Explicit copy
    Grid3D clone() const {
        Grid3D ret;
        ret.m_width  = m_width;
        ret.m_height = m_height;
        ret.m_depth  = m_depth;
        ret.m_data   = m_data;
        return ret;
    }

    bool is_valid() const {
        return m_data.size() == m_width * m_height * m_depth;
    }
    size_t data_size() const {
        return m_data.size();
    }

    T const& get(size_t x, size_t y, size_t z) const {
        return m_data[index(x, y, z)];
    }

    void set(size_t x, size_t y, size_t z, T value) {
        m_data[index(x, y, z)] = value;
    }

    T& operator()(size_t x, size_t y, size_t z) {
        return m_data[index(x, y, z)];
    }

    T const& operator()(size_t x, size_t y, size_t z) const {
        return m_data[index(x, y, z)];
    }

    // Getters for dimensions
    size_t width() const { return m_width; }
    size_t height() const { return m_height; }
    size_t depth() const { return m_depth; }

    // Iterate all values in the grid, with values
    template <class Function>
    void iterate(Function function) {
        for (int z = 0; z < m_depth; z++) {
            for (int y = 0; y < m_height; y++) {
                for (int x = 0; x < m_width; x++) {
                    T val = m_data[index(x, y, z)];
                    function(val, x, y, z);
                }
            }
        }
    }

    bool equals(const Grid3D<T>& other, float error) const {
        if (m_width != other.m_width ||
            m_height != other.m_height ||
            m_depth != other.m_depth) {
            return false;
        }
        for (int i=0; i < m_data.size(); i++) {
            auto value = std::abs(m_data[i] - other.m_data[i]);
            if (value > error) {
                return false;
            }
        }
        return true;
    }
};

template <class T>
class Grid2D {
private:
    size_t         m_width  = 0;
    size_t         m_height = 0;
    std::vector<T> m_data;

    // Convert 2D indices to a flat 1D index
    size_t index(size_t x, size_t y) const {
        if (x >= m_width || y >= m_height)
            throw std::out_of_range("Grid2D index out of range");
        return x + m_width * y;
    }

public:
    // Constructors
    Grid2D() { }
    Grid2D(size_t width, size_t height, T initial_value = T())
        : m_width(width),
          m_height(height),
          m_data(width * height, initial_value) { }

    float get(size_t x, size_t y) const { return m_data[index(x, y)]; }

    void set(size_t x, size_t y, float value) { m_data[index(x, y)] = value; }

    T& operator()(size_t x, size_t y) { return m_data[index(x, y)]; }

    T const& operator()(size_t x, size_t y) const {
        return m_data[index(x, y)];
    }

    // Getters for dimensions
    size_t width() const { return m_width; }
    size_t height() const { return m_height; }
};

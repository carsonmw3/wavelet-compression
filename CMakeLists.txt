cmake_minimum_required(VERSION 3.16)

project(wavelet-compression LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wavelet-compression main.cpp)

include(GNUInstallDirs)
install(TARGETS wavelet-compression
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

set(CPM_DOWNLOAD_VERSION 0.27.2)
set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")

if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake")
    file(DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake ${CPM_DOWNLOAD_LOCATION})
endif()

include(${CPM_DOWNLOAD_LOCATION})

find_package(PkgConfig REQUIRED)
pkg_check_modules(LZMA REQUIRED liblzma)

# include_directories(${LZMA_INCLUDE_DIRS})
# link_directories(${LZMA_LIBRARY_DIRS})
# target_link_libraries(wavelet-compression ${LZMA_LIBRARIES})

# This adds include directories and compiler flags
target_include_directories(wavelet-compression PRIVATE ${LZMA_INCLUDE_DIRS})
target_compile_options(wavelet-compression PRIVATE ${LZMA_CFLAGS_OTHER})

# This adds the linker flags
target_link_directories(wavelet-compression PRIVATE ${LZMA_LIBRARY_DIRS})
target_link_libraries(wavelet-compression PRIVATE ${LZMA_LIBRARIES})

# Link filesystem library explicitly if needed (GCC <= 8)
target_link_libraries(wavelet-compression PRIVATE stdc++fs)

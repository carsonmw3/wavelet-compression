cmake_minimum_required(VERSION 3.16)

set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/third_party/)

list(PREPEND CMAKE_PREFIX_PATH ${THIRD_PARTY_DIR})

project(wavelet-compression LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wavelet-compression main.cpp
    grid.h
    preprocess.cpp
    preprocess.h
    readandwrite.cpp
    readandwrite.h
    box-structs.h
    compressor.cpp
    compressor.h
    decompressor.cpp
    decompressor.h
    calc-loss.cpp
    calc-loss.h)

include(GNUInstallDirs)
install(TARGETS wavelet-compression
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


find_package(LibLZMA REQUIRED)
target_link_libraries(wavelet-compression PRIVATE LibLZMA::LibLZMA)

find_package(AMReX REQUIRED)
target_link_libraries(wavelet-compression PRIVATE AMReX::amrex_3d)

# Link filesystem library explicitly if needed (GCC <= 8)
if (CMAKE_COMPILER_IS_GNUCXX)
    target_link_libraries(wavelet-compression PRIVATE stdc++fs)
endif()

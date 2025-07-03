cmake_minimum_required(VERSION 3.16)

set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/third_party/)

list(PREPEND CMAKE_PREFIX_PATH ${THIRD_PARTY_DIR})

project(wavelet-compression LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    src/grid.h
    src/preprocess.cpp
    src/preprocess.h
    src/readandwrite.cpp
    src/readandwrite.h
    src/box-structs.h
    src/compressor.cpp
    src/compressor.h
    src/decompressor.cpp
    src/decompressor.h
    src/calc-loss.cpp
    src/calc-loss.h
    src/writeplotfile.cpp
    src/writeplotfile.h
    src/iterator.h
    src/argparse.cpp
    src/argparse.h
    src/tmpdir.h
    src/modes.h
    src/modes.cpp
)

add_executable(wavelet-compression src/main.cpp ${SOURCES})

target_compile_definitions(wavelet-compression PRIVATE DOCTEST_CONFIG_DISABLE)

add_executable(wavelet-compression-tests src/test.cpp ${SOURCES})

include(GNUInstallDirs)
install(TARGETS wavelet-compression
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


find_package(LibLZMA REQUIRED)
target_link_libraries(wavelet-compression PRIVATE LibLZMA::LibLZMA)

find_package(AMReX REQUIRED)
target_link_libraries(wavelet-compression PRIVATE AMReX::amrex_3d)

find_package(doctest REQUIRED)
target_link_libraries(wavelet-compression-tests PRIVATE
    AMReX::amrex_3d
    LibLZMA::LibLZMA
    doctest::doctest
)

# Link filesystem library explicitly if needed (GCC <= 8)
if (CMAKE_COMPILER_IS_GNUCXX)
    target_link_libraries(wavelet-compression PRIVATE stdc++fs)
    target_link_libraries(wavelet-compression-tests PRIVATE stdc++fs)
endif()

if  (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(wavelet-compression PUBLIC -fsanitize=address)
    target_link_options(wavelet-compression PUBLIC -fsanitize=address)
    target_compile_options(wavelet-compression-tests PUBLIC -fsanitize=address)
    target_link_options(wavelet-compression-tests PUBLIC -fsanitize=address)
endif()

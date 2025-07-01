# Wavelet-based compression for AMR data
A tool to compress Adaptive Mesh Refinement data in hopes of being able to store more of it. Allows for compression/loss tradeoffs depending on needs.

## Usage
There are three main modes that `./wavelet-compression` can be run from:
- Compression mode with `-c`: Compresses data into a specified directory
- Decompression mode (no flag): Decompresses data compressed with `./wavelet-compression` from a specified directory
- Estimate mode with `-estimate`: runs a short test with a limited amount of data to give an estimate of the following compression metrics:
  - (1) Root Mean Squared Error (RMSE) of original data vs. decompressed data, i.e. on average, the difference of any given data value across the dataset from its corresponding original after compression/decompression,
  - (2) Adjusted loss metric (RMSE divided by data range),
  - (3) compressed size (as percentage of original size).

For the compression and estimate modes, prior to the flag, the user must specify:
- `datadir`, the path to the raw AMReX plotfiles to be compressed
- `mintime`, the lowest (inclusive) timestep being compressed
- `maxtime`, the highest (inclusive) timestep being compressed
- `minlevel`, the lowest (inclusive) refinement level being compressed
- `maxlevel`, the highest (inclusive) refinement level being compressed
- `components` to be compressed; look in AMReX Header files to determine the indices of desired component(s)
- `keep`, the percentage used to calculate the threshold for keeping wavelet coefficients. Higher values lead to less compression but higher accuracy/lower loss. To start, try keep=0.99, 0.999, 0.9999 in `-estimate` mode.
- `compressedDir`, the directory where the compressed data will be written (must already exist)

For the decompression mode, the user need only specify `compressedDir`.

Thus, an example run could look like:
`./wavelet-compression datadir="../../../combustiondata/" mintime=74 maxtime=79 minlevel=0 maxlevel=3 components=0 6 25 46 keep=0.999 compressedDir="../../wavelet/" -c` 

or `./wavelet-compression datadir="../../../combustiondata/" mintime=74 maxtime=74 minlevel=0 maxlevel=1 components=6 keep=0.9999 compressedDir="../../wavelet/" -estimate`

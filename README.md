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
- `mintime`, the name of the directory in `datadir` containing the first (inclusive) timestep being compressed
- `maxtime`, the name of the directory in `datadir` containing highest (inclusive) timestep being compressed
  - Note: `mintime` and `maxtime` should be the full name of the directory, including non-digit characters. The code will automatically look for analagous directories with names containing numbers that fall between `mintime` and `maxtime`, and indicate which files will be compressed. For example, if `mintime="plt07400"` and `maxtime="plt07800"`, any directory with a similarly formatted name ending with a number between 7400 and 7800 will be compressed.
- `minlevel`, the lowest (inclusive) refinement level being compressed
- `maxlevel`, the highest (inclusive) refinement level being compressed
- `components` to be compressed, in the format `"{component_name} {component_name} ..."`; look in AMReX Header files to determine the name of components, as input must match names in Header files exactly
- `keep`, the percentage used to calculate the threshold for keeping wavelet coefficients. Higher values lead to less compression but higher accuracy/lower loss. To start, try keep=0.99, 0.999, 0.9999 in `-estimate` mode.
- `compressedDir`, the directory where the compressed data will be written (must already exist)

For the decompression mode, the user need only specify `compressedDir`.

Thus, an example run could look like:
`./wavelet-compression datadir="../../../combustiondata/" mintime=74 maxtime=79 minlevel=0 maxlevel=3 components="density Temp pressure x_velocity" keep=0.999 compressedDir="../../wavelet/" -c` 

or `./wavelet-compression datadir="../../../combustiondata/" mintime=74 maxtime=74 minlevel=0 maxlevel=1 components="Y(CH2O)" keep=0.9999 compressedDir="../../wavelet/" -estimate`

# Wavelet-based compression for AMR data
A tool to compress Adaptive Mesh Refinement data in hopes of being able to store more of it. Allows for compression/loss tradeoffs depending on needs.

## Installation
Downloading and compiling the required packages for this tool can be done by running the `install_deps.py` script, which will create a `third_party` directory.

After that, to build `wavelet-compression`:
- `mkdir build`
- `cd build`
- `cmake <path_to_wavelet-compression>`
- `make`

You should then be able to run `./wavelet-compression` followed by the options described below.


## Usage
There are three main modes that `./wavelet-compression` can be run from:
- Compression mode with `-c`: Compresses data into a specified directory
- Decompression mode with `-d`: Decompresses data compressed with `./wavelet-compression` from a specified directory
- Estimate mode with `-estimate`: runs a short test with a limited amount of data to give an estimate of the following compression metrics:
  - (1) Root Mean Squared Error (RMSE) of original data vs. decompressed data, i.e. on average, the difference of any given data value across the dataset from its corresponding original after compression/decompression,
  - (2) Adjusted loss metric (RMSE divided by data range),
  - (3) compressed size (as percentage of original size).

For the compression and estimate modes, prior to the flag, the user must specify:
- `datadir`, the path to the raw AMReX plotfiles to be compressed
- `minfile`, the name of the directory in `datadir` containing the first (inclusive) timestep being compressed
- `maxfile`, the name of the directory in `datadir` containing highest (inclusive) timestep being compressed
  - Note: `minfile` and `maxfile` should be the full name of the directory, including non-digit characters. The code will automatically look for analagous directories with names containing numbers that fall between `minfile` and `maxfile`, and indicate which files will be compressed. For example, if `minfile="plt07400"` and `maxfile="plt07800"`, any directory with a similarly formatted name ending with a number between 7400 and 7800 will be compressed.
- `minlevel`, the lowest (inclusive) refinement level being compressed
- `maxlevel`, the highest (inclusive) refinement level being compressed
- `components` to be compressed, in the format `"{component_name} {component_name} ..."`; look in AMReX Header files (e.g. `plt07400/Header`) to determine the name of components, as input must match names in Header files exactly. Also, components must be listed in the same order they appear in the Header file.
- `keep`, the percentage used to calculate the threshold for keeping wavelet coefficients. Higher values lead to less compression but higher accuracy/lower loss. To start, try keep=0.99, 0.999, 0.9999 in `-estimate` mode.
- `compressedDir`, the directory where the compressed data will be written

For the decompression mode, the user need only specify `compressedDir`.

Thus, an example run could look like:
`./wavelet-compression datadir="../../../combustiondata/" minfile="plt07400" maxfile="plt07900" minlevel=0 maxlevel=3 components="density Temp pressure x_velocity" keep=0.999 compressedDir="../../wavelet/" -c` 

or `./wavelet-compression datadir="../../../combustiondata/" minfile="plt07400" maxfile="plt07400" minlevel=0 maxlevel=1 components="Y(CH2O)" keep=0.9999 compressedDir="../../wavelet/" -estimate`

or `./wavelet-compression compresseddir="../../wavelet/" out="../../regenerated-plotfiles/" -d`

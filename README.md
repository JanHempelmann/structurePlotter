# structurePlotter

Just a small program that converts a .vesta file into a 2D-projected pdf structure.

## Dependencies

- Eigen library for linear algebra stuff ([eigen.tuxfamily.org/](https://eigen.tuxfamily.org/index.php?title=Main_Page))
- Boost library for pointer handling (https://www.boost.org/)
- Cairo for PDF/SVG rendering (https://www.cairographics.org/)

## Building on macOS

Without Homebrew/MacPorts, the easiest way to get the dependencies is via a conda environment (using Miniconda/Anaconda):

```sh
conda create -n structureplotter -c conda-forge cmake eigen boost cairo pkg-config zlib freetype cxx-compiler
conda activate structureplotter
cmake -S . -B build
cmake --build build
```

The resulting binary is `build/structurePlotter`. Run it with one or more `.vesta` files as arguments:

```sh
./build/structurePlotter path/to/structure.vesta
```

If you have Homebrew instead, `brew install cmake eigen boost cairo pkg-config` followed by the same `cmake` commands should also work.

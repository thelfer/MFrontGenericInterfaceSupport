Basic Installation
==================

These are generic installation instructions.

cmake attempts to guess correct values for various system-dependent
variables used during compilation.  It uses those values to create 
a `Makefile' in each directory of the package.  It may also create one
or more `.h' files containing system-dependent definitions..

The simplest way to compile this package is:

1. `cd` to the directory containing the package's source code and type
  `cmake .` to configure the package for your system. This will build
  `MGIS` in place.

   However, we strongly recommend building `MGIS` out of the source
   tree as follows:
   
   - Create a `build` directory et go in this directory
   - Call `cmake <path_to_MGIS_SOURCES>` to configure the directory

   Useful options are listed below.

   Running `cmake` takes awhile.  While running, it prints some
   messages telling which features it is checking for.

2. Type `make' to compile the package.

3. Optionally, type `make check' to run any self-tests that come with
  the package.

4. Type `make install' to install the programs and any data files and
  documentation.

5. You can remove the program binaries and object files from the source
  code directory by typing `make clean'.

Options
=======

- `enable-c-bindings`: compiles the `C` bindings (default=OFF)
- `enable-fortran-bindings`: compiles bindings for the `Fortran2003`
  language (default=OFF)
- `enable-python-bindings`: compiles the `Python` bindings. This requires
  the `Boost/Python` library to be available (default=OFF)
- `enable-fenics-bindings`: compiles the `FEniCS` bindings. Those
  bindings are experimental and very limited. To use `MGIS` with
  `FEniCS`, usage of the `Python` bindings are encouraged.
- `enable-julia-bindings`: compiles the `Julia` bindings. This requires
  the `CxxWrap` library to be available.
- `enable-website`: generate the TFEL website (enabled by default if
  pandoc is found)
- `enable-portable-build`: do not use processor specific flags.
- `enable-static`: compiles static libraries
- `enable-doxygen-doc`: enable the generation of the API documentation
  using with `Doxygen`.

cmake usefull variables
=======================

- `CMAKE_BUILD_TYPE`           : two values are supported 'Release' and 'Debug'
- `CASTEM_INSTALL_PATH`        : specify where the castem has been installed
- `CMAKE_TOOLCHAIN_FILE`       : specify a tool chain file
- `Python_ADDITIONAL_VERSIONS` : select the `python` version to
  			         use. Note that only the major and minor
  			         version of python shall be passed, not
  			         the revision version or the detection
  			         fails.

cmake typical usage
===================

mkdir build
cd build
cmake ../MFrontGenericInterfaceSupport/ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/th202608/codes/MFrontGenericInterfaceSupport/master/install-python-3.5/ -Denable-c-bindings=ON -Denable-fortran-bindings=ON  -Denable-python-bindings=ON -DPython_ADDITIONAL_VERSIONS=3.5
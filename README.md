# TIOGA - Topology Independent Overset Grid Assembler

[![TIOGA CI](https://github.com/Exawind/tioga/actions/workflows/ci.yml/badge.svg)](https://github.com/Exawind/tioga/actions/workflows/ci.yml)

TIOGA is a library for overset grid assembly on parallel distributed systems
Copyright (C) 2015 Jay Sitaraman

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Contact: Jay Sitaraman, jsitaraman@gmail.com, (301) 741 3216, Parallel Geometric Algorithms LLC, 814 S Mary Ave, Sunnyvale, CA 94087

## Functionality

TIOGA can perform overset grid connectivity in 3-D between
multiple unstructured (or structured) meshes that are in a distributed
computing environment, i.e. each mesh is partitioned in to multiple
parts. It can accept high-order call-back functions to perform p-consistent
interpolation and searches for formulations with internal degrees of freedom
within a computational element.

## Notes

TIOGA is free software since it was developed in the personal
time of the author. It is expected to serve as an academic/research
counterpart for PUNDIT (which is the product of the CREATE A/V program
and is export controlled). TIOGA has a subset of the functionality of
PUNDIT and is about 3x slower owing to the use of Alternating Digital Tree (ADT)
searches as the baseline point-location algorithm.

## News

TIOGA is currently under development with funding from the DoE ExaWind
program towards developing overset capability in the DoE ExaWind codes
for large scale wind farm simulations. The original code is available
[here](https://github.com/jsitaraman/tioga).

## References

Michael J. Brazell, Jayanarayanan Sitaraman, Dimitri J. Mavriplis, An overset mesh approach for 3D mixed element
high-order discretizations, In Journal of Computational Physics, Volume 322, 2016,
Pages 33-51, ISSN 0021-9991, https://doi.org/10.1016/j.jcp.2016.06.031.
(http://www.sciencedirect.com/science/article/pii/S002199911630256X)


Roget, B. and Sitaraman, J., "Robust and efficient overset grid assembly for partitioned unstructured meshes",
Journal of Computational Physics, v 260, March 2014, Pages 1-24

Brazell, M., Sitaraman J. and Mavriplis D.,"An Overset Mesh Approach for 3D Mixed
Element High Order Discretizations", Proceedings of 2014 Overset Grid Symposium,
Atlanta, GA, Oct 6-9, 2014.
http://2014.oversetgridsymposium.org/assets/presentations/3_1/Brazell_ogs_2014.pdf

## Building

### Building and installing TIOGA using CMake

TIOGA has been configured to use CMake to configure, build, and install the
library. The user can choose standard CMake options as well as additional
TIOGA-specific options to customize the build and installation process. A brief
description of the CMake-based build process and configuration options are
described in this document.

The minimal dependencies to build TIOGA on your system are CMake, a working C,
C++, and Fortran compilers as well as an MPI library along with its headers. If
the dependencies are satisfied, then execute the following commands to clone and
build TIOGA:

```
git clone <TIOGA_GITHUB_URL>
cd tioga

# Create a build directory
mkdir build

# Configure build using auto-discovered parameters
cmake ../

# Build the library
make
```

When the steps are successfully executed, the compiled static library is located
in `tioga/build/src/libtioga.a`.

#### Building `driver` and `gridGen` executables

By default, CMake does not build the `tioga.exe` driver code or the `buildGrid`
executable. To enable these at configure phase:

```
cmake -DBUILD_TIOGA_EXE:BOOL=ON -DBUILD_GRIDGEN_EXE:BOOL=ON ../
```

followed by `make`. The executables will be located in `build/driver/tioga.exe`
and `build/gridGen/buildGrid` respectively.

#### Customizing compilers

To use different compilers other than what is detected by CMake use the
following configure command:

```
CC=mpicc CXX=mpicxx FC=mpif90 cmake ../
```

#### Release, Debug, and other compilation options

Use `-DCMAKE_BUILD_TYPE` with `Release`, `Debug` or `RelWithDebInfo` to build
with different optimization or debugging flags. For example,

```
cmake -DCMAKE_BUILD_TYPE=Release ../
```

You can also use `CMAKE_CXX_FLAGS`, `CMAKE_C_FLAGS`, and `CMAKE_Fortran_FLAGS`
to specify additional compile time flags of your choosing. For example,

```
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_Fortran_FLAGS="-fbounds-check -fbacktrace" \
  ../
```

##### Custom install location

Finally, it is usually desirable to specify the install location when using
`make install` when using TIOGA with other codes.

```
# Configure TIOGA several options
CC=mpicc CXX=mpicxx FC=mpif90 cmake \
  -DCMAKE_INSTALL_PREFIX=${HOME}/software/ \
  -DBUILD_TIOGA_EXE=ON \
  -DBUILD_GRIDGEN_EXE=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_Fortran_FLAGS="-fbounds-check -fbacktrace" \
  ../

# Compile library and install at user-defined location
make && make install
```

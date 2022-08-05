`Documentation in progress...`

## Overview

Kinematic and vertex fitting is widely used data analysis technique in particle physics experiments. This technique can be used in order to separate events corresponding to different kinematic hypotheses and to restore interaction and decay vertices. In addition, the use of this technique often leads to significant improvements in resolution of some quantities, such as various invariant masses.

`KFBase` contains the implementation of base classes for kinematic and vertex fitting. This package was developed for the CMD-3 experiment. However, the `KFBase` package is not detector-dependent and can be used in another experiment. The detector-dependent part of the kinematic and vertex fitting software for the CMD-3 experiment is contained in package `KFCmd`: https://github.com/sergeigribanov/KFCmd.

## Dependencies
1. `Eigen 3` template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms: https://eigen.tuxfamily.org/index.php?title=Main_Page (successfully tested with version `3.4`).
2. `ROOT` data analysis framework used by high energy physics and others: https://root.cern (successfully tested with version `6.26`).

## Installation
1. Get the source code:
```console
git clone https://github.com/sergeigribanov/KFBase
```
2. Setup `ROOT` environment:
```console
source <path to ROOT installation>/bin/thisroot.sh
```
3. Create a build directory:
```console
mkdir <path to a build directory>
cd <path to a build directory>
```
4. Run CMake:
```console
cmake -DCMAKE_INSTALL_PREFIX=<installation prefix> <path to the source code directory>
```
To set `C++` standard, `-DCMAKE_CXX_STANDARD` option can be used:
```console
cmake -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=<installation prefix> <path to the source code directory>
``` 
The `KFBase` package has been successfully tested with both `C++11` and `C++17` standards. The same standard should be used that was used to build the `ROOT` framework.

In some cases, the option `-DEigen3_DIR` may be required in order to find `Eigen 3` installation. See page https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html for more details.

5. Build the package:
```console
make
```
6. Install the package:
```console
make install
```

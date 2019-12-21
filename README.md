# KFBase

A library contained abstract classes for kinemtaic fit.

## Installing

1. git clone https://github.com/sergeigribanov/KFBase.git
2. Create a directory to build the package in a suitable location and change the current directory to this one.
3. To build the package run the following commands: 
    1. cmake -DCMAKE_INSTALL_PREFIX=\<KFBase installation prefix\> -DEIGEN3_INCLUDE_DIR=\<path to Eigen3 installation\> -DCCGO_DIR=\<path to CCGO installation\> \<path to KFBase source code\>
    2. make
    3. make install

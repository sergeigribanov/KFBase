# KFBase

A library that contains abstract classes for kinemtaic fit.

## Dependencies

1. Eigen3: http://eigen.tuxfamily.org/index.php?title=Main_Page
2. ROOT: https://root.cern.ch
3. CCGO: https://github.com/sergeigribanov/ccgo

## Installing

1. git clone https://github.com/sergeigribanov/KFBase.git
2. Create a directory to build the package in a suitable location and change the current directory to this one.
3. Setup ROOT environment.
4. Make sure that CCGO package has been installed.
5. To build the package run the following commands: 
    1. cmake -DCMAKE_INSTALL_PREFIX=\<KFBase installation prefix\> \<path to KFBase source code\>
    2. make
    3. make install

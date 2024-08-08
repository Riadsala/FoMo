#!/bin/bash
# File: spatialWrapper.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -l mem_free=10G
#$ -m be
#$ -M ah19679@essex.ac.uk
#$ -o output_spatial.txt
export PATH=/usr/local/gcc9.3/bin/:$PATH
export LD_LIBRARY_PATH=/usr/local/gcc9.3/lib64:$LD_LIBRARY_PATH
export CC=/usr/local/gcc9.3/bin/gcc
export CPP=/usr/local/gcc9.3/bin/cpp
export CXX=/usr/local/gcc9.3/bin/c++
CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS} -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H"
export PATH=/usr/local/gdal3/bin/:/usr/local/poppler-0.47.0/bin/:/usr/local/ncdf4.8.0/bin/:/usr/local/hdf5-1.12.0/bin/:$PATH
export LD_LIBRARY_PATH=/usr/local/gdal3/lib/:/usr/local/poppler-0.47.0/lib/:/usr/local/ncdf4.8.0/lib/:/usr/local/mpfr-4.0.1/lib/:/usr/local/hdf5-1.12.0/lib/:$LD_LIBRARY_PATH:/usr/local/proj6/lib/:/usr/lib64/

./cluster_model_fitting_test.R

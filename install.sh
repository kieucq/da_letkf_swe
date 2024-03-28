#!/bin/sh
# NOTE: This sciprt compiles the whle SWE-LETKF system. It requires 
#       a fortran compiler and several internal libs. The only Fortran
#       compilers that havve been tested so far are ifort and PGI.
#
# AUTH: Chanh Kieu (ckieu@iu.edu)
#=====================================================================
set -x
MAIN_DIR=/N/u/ckieu/Karst/model/da_letkf_swe
#FC="ifort"
#MFC="/N/u/ckieu/Quartz/local/libs/openmpi-4.0.1-intel/bin/mpif90"
FC="pgf90"
MFC="mpif90"
FCFLAG="-traceback"
INC="$MAIN_DIR/registry"
MPI="No"
if [ $MPI == "Yes" ] && [ $FC != "pgf90" ]; then
   echo "MPI mode can only work with PGI ...stop"
   exit 1
fi
#
# compiling diagnostic directory
#
echo "Compiling in dir: dig"
cd $MAIN_DIR/dig
rm -rf *.exe *.o *.mod
$FC $FCFLAG -o ./ana.exe ./ana.f90
#
# compiling letkf directory
#
echo "Compiling in dir: letkf"
cd $MAIN_DIR/letkf
rm -rf *.exe *.o *.mod
$FC $FCFLAG -c mt19937ar.f90
$FC $FCFLAG -c common.f90
$FC $FCFLAG -c netlib.f
$FC $FCFLAG -c matrix_inv.f90 
$FC $FCFLAG -c common_mtx.f90
if [ $MPI == "Yes" ]; then
   $MFC $FCFLAG -c common_mpi.f90
   $MFC $FCFLAG -o letkf_mpi.exe letkf_mpi.f90 common_mtx.f90      \ 
                   common.f90 mt19937ar.f90 netlib.f               \
                   matrix_inv.f90 common_mpi.f90
else
   $FC $FCFLAG -o letkf_serial.exe letkf_serial.f90 common_mtx.f90 \
                  common.f90 mt19937ar.f90 netlib.f matrix_inv.f90
fi
#
# compiling ctl directory
#
echo "Compiling in dir: ctl"
cd $MAIN_DIR/ctl
$FC $FCFLAG -o ctl.exe swe.f
#
# compiling obs directory
#
echo "Compiling in dir: obs"
cd $MAIN_DIR/obs
rm -rf *.exe *.o *.mod
$FC $FCFLAG -c mt19937ar.f90
$FC $FCFLAG -o obs.exe obs.f90 mt19937ar.f90
#
# compiling truth directory
#
echo "Compiling in dir: truth"
cd $MAIN_DIR/truth
rm -rf *.exe *.o *.mod
$FC $FCFLAG -o truth.exe swe.f 
$FC $FCFLAG -o bvortex.exe bvortex.f90
#
# compiling model directory
#
echo "Compiling in dir: model"
cd $MAIN_DIR/model
rm -rf *.exe *.o *.mod
$FC $FCFLAG -o swe.exe swe.f
#
# compiling ultility directory
#
echo "Compiling in dir: utils"
cd $MAIN_DIR/utils
rm -rf *.exe *.o *.mod
$FC $FCFLAG -o  mean.exe mean.f90
#
# compiling ini directory
#
echo "Compiling in dir: ini"
cd $MAIN_DIR/ini
rm -rf *.exe *.o *.mod
$FC $FCFLAG -c mt19937ar.f90
$FC $FCFLAG -o ini_letkf.exe ini_letkf.f90 mt19937ar.f90
$FC $FCFLAG -o bgd.exe bgd.f
cd ../
ls -la */*.exe
check_compile=`ls -la */*.exe | wc -l`
if [ $check_compile == 10 ]; then
    echo "Compilation is successful"
else
    echo "Compilation failed. Check again"
fi
echo "DONE INSTALLING"

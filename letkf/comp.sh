#mpif90 -o letkf_mpi.exe letkf_mpi.f90 common_mtx.f90 common.f90 mt19937ar.f90 netlib.f matrix_inv.f90 common_mpi.f90 module_interface.f90
#mpif90 -o 4dletkf_mpi.exe 4dletkf.f90 common_mtx.f90 common.f90 mt19937ar.f90 netlib.f matrix_inv.f90 common_mpi.f90 module_interface.f90 module_4dletkf.f90 module_utils.f90
pgf90 -o letkf_serial.exe letkf_serial.f90 common_mtx.f90 common.f90 mt19937ar.f90 netlib.f matrix_inv.f90


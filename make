#... with full debug options
#h5fc -g -traceback -mkl -warn unused mod_phdf5.f90 mod_blocks.f90 test_matrix_c.f90 -o test_matrix_c
h5fc -g -traceback -mkl -warn unused mod_phdf5.f90 mod_blocks.f90 test_read_c.f90 -o test_read_c

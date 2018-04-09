program main
  use hdf5
  use mod_phdf5
  use mod_blocks

  implicit none
  type(block1d) :: testblock
  integer :: nblocks_, blocksize_
  integer :: k
  ! set global parameters
  nblocks_=
  blocksize_=
  ! global set-up for hdf5
  call 
  ! generate a test-matrix, distribute it over blocks and write to HDF5
  do k=1, nblocks_
    ! set up testblock
    testblock%nblocks=nblocks_
    testblock%blocksize=blocksize_
    testblock%il=(k-1)*blocksize_+1
    testblock%iu=k*blocksize_
    testblock%offset=(k-1)*blocksize_
    testblock%id=k
    ! allocate and fill zcontent
    if (allocated(testblock%zcontent)) deallocate(testblock%zcontent)
    allocate(testblock%zcontent(blocksize_))
    testblock%zcontent(:)=cmplx(float(k),0.0d0)
    ! write block to hdf5 
    call put_block1d(testblock,.True.,dataset_id) 
  end do
end program


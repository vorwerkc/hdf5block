program main
  use hdf5
  use mod_phdf5
  use mod_blocks

  implicit none
  type(block1d) :: testblock
  character(1024) :: fname
  integer :: nblocks_, blocksize_, matsize_
  integer :: k
  integer(hid_t) :: file_id, dataset_id
  character(1024) :: path, groupname
  ! set global parameters
  fname='vector.h5'
  nblocks_=4 ! number of blocks
  blocksize_=10 ! blocksize
  matsize_=nblocks_*blocksize_ ! total matrix size
  path='/'
  groupname='A'
  ! global set-up for hdf5
  call phdf5_initialize(fname,file_id)
  ! matrix set-up for hdf5
  call phdf5_setup(matsize_,.False.,groupname,path,file_id,dataset_id) 
  print *, 'dataset=', dataset_id
  
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
    if (allocated(testblock%dcontent)) deallocate(testblock%dcontent)
    allocate(testblock%dcontent(blocksize_))
    testblock%dcontent(:)=float(k)
    ! write block to hdf5 
    call put_block1d(testblock,dataset_id) 
  end do
  ! finalize
  call phdf5_cleanup(dataset_id)
  call phdf5_finalize(file_id)
end program


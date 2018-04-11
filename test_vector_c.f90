program main
  use hdf5
  use mod_phdf5
  use mod_blocks

  implicit none
  type(block1d) :: testblock, readblock
  character(1024) :: fname
  integer :: nblocks_, blocksize_, matsize_, dimsg_(1)
  integer :: k, l
  integer(hid_t) :: file_id, dataset_id
  character(1024) :: path, groupname
  ! set global parameters
  fname='vector.h5'
  nblocks_=4 ! number of blocks
  blocksize_=10 ! blocksize
  matsize_=nblocks_*blocksize_ ! total matrix size
  dimsg_(1)=matsize_
  path='/'
  groupname='A'
  ! global set-up for hdf5
  call phdf5_initialize(fname,file_id)
  ! matrix set-up for hdf5
  call phdf5_setup_write(1,dimsg_,.True.,groupname,path,file_id,dataset_id) 
  
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
    call put_block1d(testblock,dataset_id) 
  end do
  ! finalize
  call phdf5_cleanup(dataset_id)

  ! open file again and read the vector to blocks
  call phdf5_setup_read(matsize_,.True.,groupname,path,file_id,dataset_id)
  do k=1,nblocks_
    ! set up the readblocks
    readblock%nblocks=nblocks_
    readblock%blocksize=blocksize_
    readblock%il=(k-1)*blocksize_+1
    readblock%iu=k*blocksize_
    readblock%offset=(k-1)*blocksize_
    readblock%id=k
    ! allocate output array for the data
    if (allocated(readblock%zcontent)) deallocate(readblock%zcontent)
    allocate(readblock%zcontent(blocksize_))
    call get_block1d(readblock,dataset_id)
    ! print the content of blocks
    print *, 'k=', k
    do l=1, blocksize_
      print *, 'readblock%zcontent(',l,')=', readblock%zcontent(l)
    end do
  end do
  ! finalize
  call phdf5_cleanup(dataset_id)
  call phdf5_finalize(file_id)
end program


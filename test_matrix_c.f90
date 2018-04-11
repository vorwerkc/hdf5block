program main
  use hdf5
  use mod_phdf5
  use mod_blocks

  implicit none
  type(block2d) :: testblock, readblock
  character(1024) :: fname
  integer :: nblocks_, blocksize_, matsize_, dimsg_(2)
  integer :: k1, k2, l1, l2
  integer(hid_t) :: file_id, dataset_id
  character(1024) :: path, groupname
  ! set global parameters
  fname='matrix.h5'
  nblocks_=2 ! number of blocks
  blocksize_=2 ! blocksize
  matsize_=nblocks_*blocksize_ ! total matrix size
  dimsg_=(/matsize_, matsize_/)
  path='/'
  groupname='A'
  ! global set-up for hdf5
  call phdf5_initialize(fname,file_id)
  ! matrix set-up for hdf5
  call phdf5_setup_write(2,dimsg_,.True.,groupname,path,file_id,dataset_id) 
  
  ! generate a test-matrix, distribute it over blocks and write to HDF5
  do k1=1, nblocks_
    do k2=1, nblocks_
      ! set up testblock
      testblock%nblocks=nblocks_
      testblock%blocksize=blocksize_
      testblock%il=(k1-1)*blocksize_+1
      testblock%iu=k1*blocksize_
      testblock%jl=(k2-1)*blocksize_+1
      testblock%ju=k2*blocksize_
      testblock%offset=(/ (k1-1)*blocksize_, (k2-1)*blocksize_ /)
      testblock%id=(/ k1, k2 /)
      ! allocate and fill zcontent
      if (allocated(testblock%zcontent)) deallocate(testblock%zcontent)
      allocate(testblock%zcontent(blocksize_, blocksize_))
      do l1=1, blocksize_
        do l2=1, blocksize_
          testblock%zcontent(l1,l2)=cmplx(float(l1),float(l2))
        end do
      end do
      ! write block to hdf5 
      call put_block2d(testblock,dataset_id)
    end do 
  end do

  ! finalize
  call phdf5_cleanup(dataset_id)

  ! open file again and read the vector to blocks
  call phdf5_setup_read(matsize_,.True.,groupname,path,file_id,dataset_id)
  do k1=1,nblocks_
    do k2=1,nblocks_
      ! set up the readblocks
      readblock%nblocks=nblocks_
      readblock%blocksize=blocksize_
      readblock%il=(k1-1)*blocksize_+1
      readblock%iu=k1*blocksize_
      readblock%jl=(k2-1)*blocksize_+1
      readblock%ju=k2*blocksize_
      readblock%offset=(/ (k1-1)*blocksize_, (k2-1)*blocksize_ /)
      readblock%id=(/ k1, k2 /)
      ! allocate output array for the data
      if (allocated(readblock%zcontent)) deallocate(readblock%zcontent)
      allocate(readblock%zcontent(blocksize_, blocksize_))
      call get_block2d(readblock,dataset_id)
      ! print the content of blocks
      print *, 'id=', readblock%id
      print *, readblock%zcontent(1,1),readblock%zcontent(1,2)
      print *, readblock%zcontent(2,1),readblock%zcontent(2,2)
    end do
  end do
  ! finalize
  call phdf5_cleanup(dataset_id)
  call phdf5_finalize(file_id)
end program


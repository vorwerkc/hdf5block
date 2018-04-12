program main
  use hdf5
  use mod_phdf5
  use mod_blocks
  use mpi

  implicit none
  type(block1d) :: testblock, readblock
  character(1024) :: fname
  integer :: nblocks_, blocksize_, matsize_, dimsg_(1)
  integer :: k, l
  integer(hid_t) :: file_id, dataset_id
  character(1024) :: path, datasetname, groupname
  complex(8), allocatable :: global(:)
  ! mpi variables
  integer :: mpierror
  integer :: comm, info
  integer :: mpi_size, mpi_rank
  integer :: start, ende, perrank
  ! set global parameters
  fname='vector.h5'
  nblocks_=4 ! number of blocks
  blocksize_=10 ! blocksize
  matsize_=nblocks_*blocksize_ ! total matrix size
  dimsg_(1)=matsize_
  path='/'
  datasetname='A'
  ! set global MPI settings
  comm=MPI_COMM_WORLD
  info=MPI_INFO_NULL
  ! MPI initialization
  call mpi_init(mpierror)
  call mpi_comm_size(comm,mpi_size,mpierror)
  call mpi_comm_rank(comm, mpi_rank,mpierror)
  ! global set-up for hdf5
  call phdf5_initialize(fname,.True.,file_id,comm)
  ! create a group
  groupname='test'
  call phdf5_create_group(file_id,path,groupname)
  path='/'//trim(groupname)//'/'
  ! matrix set-up for hdf5
  call phdf5_setup_write(1,dimsg_,.True.,datasetname,path,file_id,dataset_id) 
  ! distribute blocks over MPI ranks
  perrank=floor(float(nblocks_)/float(mpi_size))
  start=(mpi_rank)*perrank+1
  ende=(mpi_rank+1)*perrank
  ! generate a test-matrix, distribute it over blocks and write to HDF5
  do k=start, ende
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
  ! finalize
  call phdf5_cleanup(dataset_id)
  ! now we try to read it again
  call phdf5_setup_read(1,dimsg_,.True.,datasetname,path,file_id,dataset_id)
  ! allocate global array
  allocate(global(perrank*blocksize_))
  ! loop over blocks
  do k=start, ende
    ! set up readblock
    readblock%nblocks=nblocks_
    readblock%blocksize=blocksize_
    readblock%il=(k-1)*blocksize_+1
    readblock%iu=k*blocksize_
    readblock%offset=(k-1)*blocksize_
    readblock%id=k
    if (allocated(readblock%zcontent)) deallocate(readblock%zcontent)
    allocate(readblock%zcontent(blocksize_))
    call get_block1d(readblock,.True.,dataset_id)
    do l=1, blocksize_
      global(l+(readblock%id-start)*blocksize_)=readblock%zcontent(l)
    end do
  end do
  ! write out the result
  if (mpi_rank .eq. 1) then
    print *, '******** rank=', mpi_rank, '**************'
    do l=1, size(global)
      print *, global(l)
    end do 
  end if
  deallocate(global)
  call phdf5_cleanup(dataset_id)
  call phdf5_finalize(file_id)
  call mpi_finalize(mpierror)
end program


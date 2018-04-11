program main
  use hdf5
  use mod_phdf5
  use mod_blocks
  use mpi

  implicit none
  type(block1d) :: testblock
  character(1024) :: fname
  integer :: nblocks_, blocksize_, matsize_, dimsg_(1)
  integer :: k
  integer(hid_t) :: file_id, dataset_id
  character(1024) :: path, groupname
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
  groupname='A'
  ! set global MPI settings
  comm=MPI_COMM_WORLD
  info=MPI_INFO_NULL
  ! MPI initialization
  call mpi_init(mpierror)
  call mpi_comm_size(comm,mpi_size,mpierror)
  call mpi_comm_rank(comm, mpi_rank,mpierror)
  ! global set-up for hdf5
  call phdf5_initialize(fname,.True.,file_id,comm)
  ! matrix set-up for hdf5
  call phdf5_setup_write(1,dimsg_,.True.,groupname,path,file_id,dataset_id) 
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
  call phdf5_finalize(file_id)
  call mpi_finalize(mpierror)
end program


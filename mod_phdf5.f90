module mod_phdf5
  interface phdf5_write
    module procedure phdf5_write_d, &
        &            phdf5_write_z
  end interface
  interface phdf5_read
    module procedure phdf5_read_d, &
        &            phdf5_read_z
  end interface
  public phdf5_initialize
  public phdf5_finalize
  public phdf5_setup_write
  public phdf5_cleanup
contains

!-------------------------------------------------------------------------------
  subroutine phdf5_initialize(fname,fparallel,file_id,comm)
    use hdf5
    use mpi
    implicit none
    character(1024), intent(in) :: fname
    logical, intent(in) :: fparallel
    integer(hid_t), intent(out) :: file_id
    integer, optional :: comm
    ! local variables
    integer (hid_t) :: plist_id
    integer :: ierr
    integer :: info
    character*100 :: errmsg
    ! MPI File creation & Global Access
    if (fparallel) then
      call h5open_f(ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5open_f returned ",I6)')ierr
        goto 10
      endif    
      ! set mpi info object
      info=MPI_INFO_NULL
      ! create file access property list w/ parallel IO access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5pcreate_f returned ",I6)')ierr
        goto 10
      endif    
      call h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5pset_fapl_mpio_f returned ",I6)')ierr
        goto 10
      endif    
      ! create the file collectively
      call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,file_id,ierr,access_prp=plist_id)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5fcreate_f returned ",I6)')ierr
        goto 10
      endif    
      ! close property list
      call h5pclose_f(plist_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5pclose_f returned ",I6)')ierr
        goto 10
      endif    
    ! Serial Access
    else
      call h5open_f(ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5open_f returned ",I6)')ierr
        goto 10
      endif    
      call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_initialize): h5fcreate_f returned ",I6)')ierr
        goto 10
      endif    
    end if
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  fname : ",A)')trim(fname)
    stop
  end subroutine

!-------------------------------------------------------------------------------
  subroutine phdf5_finalize(file_id)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: file_id
    ! local variable
    integer :: ierr
    character*100 :: errmsg

    call h5fclose_f(file_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_finalize): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  file_id : ",I4)')file_id
    stop
    call h5close_f(ierr)
  
  end subroutine 

!-------------------------------------------------------------------------------
  subroutine phdf5_setup_write(ndims,dims,fcomplex,dname,path,file_id,dataset_id)
    use hdf5
    implicit none
    integer, intent(in) :: ndims
    integer, dimension(ndims), intent(in) :: dims
    logical, intent(in) :: fcomplex
    character(1024), intent(in) :: dname, path
    integer(hid_t), intent(in) :: file_id
    integer(hid_t), intent(out) :: dataset_id
    !local variables
    integer :: ndims_
    ! HDF5 variables
    integer(hid_t) :: dataspace_id, group_id
    integer(hsize_t), allocatable :: dims_(:)
    integer :: ierr 
    character*100 :: errmsg
    ! if the dataset represents complex data, create 2D array
    if (fcomplex) then
      ndims_=ndims+1
      allocate(dims_(ndims_))
      dims_(1)=2
      dims_(2:)=dims
    else
      ndims_=ndims
      allocate(dims_(ndims_))
      dims_(:)=dims
    end if
    ! create the dataspace
    call h5screate_simple_f(ndims_,dims_,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_setup): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif    
    ! open group
    call h5gopen_f(file_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_setup): h5gopen_f returned ",I6)')ierr
      goto 10
    endif    
    ! create the dataset
    call h5dcreate_f(group_id,trim(dname),H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): h5dcreate_f returned ",i6)')ierr
      goto 10
    endif    
    ! close dataset, group, dataspace
    call h5gclose_f(group_id,ierr)
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): closing returned ",i6)')ierr
      goto 10
    endif    
    deallocate(dims_)
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims_
    write(*,'("  fname : ",A)')trim(dname)
    write(*,'("  path  : ",A)')trim(path)
    stop
  end subroutine
  
!-------------------------------------------------------------------------------
  subroutine phdf5_setup_read(dims,fcomplex,dname,path,file_id,dataset_id)
    use hdf5
    implicit none
    integer, intent(in) :: dims
    logical, intent(in) :: fcomplex
    character(1024), intent(in) :: dname, path
    integer(hid_t), intent(in) :: file_id
    integer(hid_t), intent(out) :: dataset_id
    !local variables
    integer :: ndims_
    ! HDF5 variables
    integer(hid_t) ::group_id
    integer(hsize_t), allocatable :: dims_(:)
    integer :: ierr 
    character*100 :: errmsg
    ! if the dataset represents complex data, create 2D array
    if (fcomplex) then
      ndims_=2
      allocate(dims_(ndims_))
      dims_=(/2, dims/)
    else
      ndims_=1
      allocate(dims_(ndims_))
      dims_(1)=dims
    end if
    ! open group
    call h5gopen_f(file_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_setup): h5gopen_f returned ",I6)')ierr
      goto 10
    endif    
    ! open the dataset
    call h5dopen_f(group_id,trim(dname),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): h5dopen_f returned ",i6)')ierr
      goto 10
    endif    
    ! close dataset, group
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): closing returned ",i6)')ierr
      goto 10
    endif    
    deallocate(dims_)
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims_
    write(*,'("  fname : ",A)')trim(dname)
    write(*,'("  path  : ",A)')trim(path)
    stop
  end subroutine
!-------------------------------------------------------------------------------
  subroutine phdf5_cleanup(dataset_id)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: dataset_id
    ! local variable
    integer :: ierr

    call h5dclose_f(dataset_id,ierr)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_write_d(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    real(8), intent(in) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(:)=dims(:)
    dimsg_(:)=dimsg(:)
    offset_(:)=offset(:)
    ! write to hdf5
    call phdf5_write_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_write_z(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    complex(8), intent(in) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)+1
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(1)=2
    dims_(2:)=dims(:)
    dimsg_(1)=2
    dimsg_(2:)=dimsg(:)
    offset_(1)=0
    offset_(2:)=offset(:)
    ! write to hdf5
    call phdf5_write_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_read_d(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    real(8), intent(out) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(:)=dims(:)
    dimsg_(:)=dimsg(:)
    offset_(:)=offset(:)
    ! write to hdf5
    call phdf5_read_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_read_z(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    complex(8), intent(out) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)+1
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(1)=2
    dims_(2:)=dims(:)
    dimsg_(1)=2
    dimsg_(2:)=dimsg(:)
    offset_(1)=0
    offset_(2:)=offset(:)
    ! write to hdf5
    call phdf5_read_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine
end module
!-----------------------------------------------------------------------------
subroutine phdf5_write_array_d(a,fparallel,ndims,dims,dimsg,offset,dataset_id)
  use hdf5
  implicit none
  real(8), intent(in) :: a(*)
  logical, intent(in) :: fparallel
  integer, intent(in) :: ndims
  integer(hsize_t), intent(in) :: dims(ndims), dimsg(ndims), offset(ndims)
  integer(hid_t), intent(in) :: dataset_id
  ! local variables
  integer :: ierr
  integer(hid_t) :: memspace_id, dataspace_id, plist_id
  ! create memoryspace
  call h5screate_simple_f(ndims,dims,memspace_id,ierr)
  ! select hyperslab in file
  call h5dget_space_f(dataset_id,dataspace_id,ierr)
  call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset,dims,ierr)
  ! write into hyperslab
  ! MPI I/O independently
  if (fparallel) then
    ! create property list for individual dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F,ierr)
    ! write the dataset
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id,plist_id)
    ! close the property list
    call h5pclose_f(plist_id,ierr)
  ! serial
  else
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id)
  end if
  ! close memory space and dataspace
  call h5sclose_f(dataspace_id,ierr)
  call h5sclose_f(memspace_id,ierr) 
end subroutine

!-----------------------------------------------------------------------------
subroutine phdf5_read_array_d(a,fparallel,ndims,dims,dimsg,offset,dataset_id)
  use hdf5
  implicit none
  real(8), intent(out) :: a(*)
  logical, intent(in) :: fparallel
  integer, intent(in) :: ndims
  integer(hsize_t), intent(in) :: dims(ndims), dimsg(ndims), offset(ndims)
  integer(hid_t), intent(in) :: dataset_id
  ! local variables
  integer :: ierr
  integer(hid_t) :: memspace_id, dataspace_id, plist_id
  ! create memoryspace
  call h5screate_simple_f(ndims,dims,memspace_id,ierr)
  ! select hyperslab in file
  call h5dget_space_f(dataset_id,dataspace_id,ierr)
  call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset,dims,ierr)
  ! read from hyperslab
  ! MPI read
  if (fparallel) then
    ! create property list for individual dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F,ierr)
    ! read dataset into memory
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id,plist_id)
    ! close the property list
    call h5pclose_f(plist_id,ierr)
  ! serial read
  else
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id)
  end if
  ! close memory space and dataspace
  call h5sclose_f(dataspace_id,ierr)
  call h5sclose_f(memspace_id,ierr) 
end subroutine

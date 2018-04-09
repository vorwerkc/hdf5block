module mod_phdf5

  public prepare_output_global
  public finalize_output_global
  public prepare_output1d
  public finalize_ouput1d
  public put_block1d
  public get_block1d
contains

!-------------------------------------------------------------------------------
  subroutine prepare_output_global(fname,file_id)
    use hdf5
    implicit none
    character(1024), intent(in) :: fname
    integer(hid_t), intent(out) :: file_id
    ! local variables
    integer :: ierr

    call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
  end subroutine

!-------------------------------------------------------------------------------
  subroutine finalize_output_global(file_id)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: file_id
    ! local variable
    integer :: ierr

    call h5fclose(file_id,ierr)
  
  end subroutine 

!-------------------------------------------------------------------------------
  subroutine prepare_output1d(dims,fcomplex,dname,path,file_id,dataset_id)
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
    integer(hid_t) :: dataspace_id, group_id
    integer(hsize_t), allocatable :: dims_(:)
    integer :: ierr 
    character*100 :: errmsg
    ! if the dataset represents complex data, create 2D array
    if (fcomplex) then
      ndims_=2
      allocate(dims_(ndims))
      dims_=(/2, dims/)
    else
      ndims_=1
      allocate(dims_(ndims))
      dims_(1)=dims
    end if
    ! create the dataspace
    call h5screate_simple_f(ndims_,dims_,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(prepare_output): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif    
    ! open group
    call h5gopen_f(file_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(prepare_output): h5gopen_f returned ",I6)')ierr
      goto 10
    endif    
    ! create the dataset
    call h5dcreate_f(group_id,trim(dname),H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(prepare_output): h5dcreate_f returned ",i6)')ierr
      goto 10
    endif    
    ! close dataset, group, dataspace
    call h5gclose_f(group_id,ierr)
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(prepare_output): closing returned ",i6)')ierr
      goto 10
    endif    
    deallocate(dims_)
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  fname : ",A)')trim(dname)
    write(*,'("  path  : ",A)')trim(path)
    stop
  end subroutine
  
!-------------------------------------------------------------------------------
  subroutine finalize_output1d(dataset_id)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: dataset_id
    ! local variable
    integer :: ierr

    call h5dclose_f(dataset_id,ierr)

  end subroutine

!-----------------------------------------------------------------------------
  subroutine put_block1d(in1d,fcomplex,dataset_id)
    use hdf5
    implicit none
    type(block1d), intent(in) :: in1d
    logical, intent(in) :: fcomplex
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer :: ndims_, ierr
    integer(hsize_t), allocatable, dimension(:) :: dims_, offset_, dimsg_
    integer(hid_t) :: memspace_id, dataspace_id
    if (fcomplex) then
      ! set dimensions of the block, since the data is complex valued, there are 2 dims
      allocate(dims_(2),offset_(2),dimsg_(2))
      ndims_=2
      dims_=(/2, in1d%blocksize/)
      ! calculate global size of array
      dimsg_=(/2, in1d%blocksize*in1d%nblocks /)
      ! set parameters of hyperslab
      offset_=(/0, in1d%offset/)
    else
      allocate(dims_(1),offset_(1),dimsg_(1))
      ndims_=1
      dims_=(/ in1d%blocksize/)
      dimsg_=(/in1d%blocksize*in1d%nblocks/)
      offset_=(/in1d%offset/)
    end if
    ! create memoryspace
    call h5screate_simple_f(ndims_,dims_,memspace_id,ierr)
    ! select hyperslab in file
    call h5dget_space_f(dataset_id,dataspace_id,ierr)
    call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset_,dims_,ierr)
    ! write into hyperslab
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,in1d%zcontent,dimsg_,ierr,file_space_id=dataspace_id,mem_space_id=memspace_id)
    ! close memory space and dataspace
    call h5sclose(dataspace_id,ierr)
    call h5sclose(memspace_id,ierr)
    ! deallocate arays
    deallocate(dims_,dimsg_,offset_)
  end subroutine
!-----------------------------------------------------------------------------
  subroutine get_block1d(in1d,fcomplex,dataset_id)
    use hdf5
    implicit none
    type(block1d), intent(inout) :: in1d
    logical, intent(in) :: fcomplex
    integer(hid_t) :: dataset_id
    ! local variables
    integer :: ndims_, ierr
    integer(hsize_t), allocatable, dimension(:) :: dims_, offset_, dimsg_
    integer(hid_t) :: memspace_id, dataspace_id
    if (fcomplex) then
      ! set dimensions of the block, since the data is complex valued, there are 2 dims
      allocate(dims_(2),dimsg_(2), offset_(2))
      ndims_=2
      dims_=(/2, in1d%blocksize/)
      ! calculate global size of array
      dimsg_=(/2, in1d%blocksize*in1d%nblocks /)
      ! set parameters of hyperslab
      offset_=(/0, in1d%offset/)
    else
      allocate(dims_(1),dimsg_(1),offset_(1))
      ndims_=1
      dims_=(/ in1d%blocksize/)
      ! calculate global size of array
      dimsg_=(/in1d%blocksize*in1d%nblocks /)
      ! set parameters of hyperslab
      offset_=(/in1d%offset/)
    end if
    ! create memoryspace
    call h5screate_simple_f(ndims_,dims_,memspace_id,ierr)
    ! select hyperslab in file
    call h5dget_space_f(dataset_id,dataspace_id,ierr)
    call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset_,dims_,ierr)
    ! allocate output
    if (allocated(in1d%zcontent)) deallocate(in1d%zcontent)
    allocate(in1d%zcontent(in1d%blocksize))
    ! read from hyperslab
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,in1d%zcontent,dimsg_,ierr,file_space_id=dataspace_id,mem_space_id=memspace_id)
    ! close memory space and dataspace
    call h5sclose(dataspace_id,ierr)
    call h5sclose(memspace_id,ierr)
    deallocate(dims_,dimsg_,offset_)
  end subroutine

end module

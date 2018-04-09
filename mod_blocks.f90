module mod_blocks
  implicit none

  type :: block1d
    integer(4) :: nblocks     ! number of blocks
    integer(4) :: blocksize   ! length of each block
    integer(4) :: nk          ! number of k-pts per block
    integer(4) :: il, iu      ! absolute first & last transition index within block
    integer(4) :: kl, ku      ! absolute first & last k-index within block
    integer(4) :: offset      ! offset
    integer(4) :: id  ! ID (in this case number) of the block
    complex(8), allocatable :: zcontent(:)
    real(8), allocatable :: dcontent(:)
  end type block1d
  
  type :: block2d
    integer(4) :: nblocks     ! number of block
    integer(4) :: blocksize   ! size of each block is blocksize x blocksize
    integer(4) :: nk          ! number of k-pts per block
    integer(4) :: il, iu, jl, ju ! absolute transition index ranges within the block  
    integer(4) :: k1l, k1u, k2l, k2u ! absolute k-index within the block
    integer(4), dimension(2) :: offset    ! offset
    integer(4), dimension(2) :: id    ! 2-D ID of block
    complex(8), allocatable :: zcontent(:,:)
    real(8), allocatable :: dcontent(:,:)
  end type block2d

  type :: block3d
    integer(4) :: nblocks                 ! number of blocks
    integer(4) :: nk                      ! number of k-pts / block
    integer(4), dimension(3) :: blocksize ! non-square 3D blocksize
    integer(4) :: kl, ku                  ! absolute k-index within block
    integer(4) :: id                      ! store the id of the block it was generated from
    complex(8), allocatable :: zcontent(:,:,:)
  end type block3d
  
  contains
!-----------------------------------------------------------------------------
  subroutine put_block1d(in1d,dataset_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block1d), intent(inout) :: in1d
    integer(hid_t) :: dataset_id
    ! local variables
    integer, dimension(1) :: dims_, dimsg_, offset_
    ! set dimension & offset
    dims_(1)=in1d%blocksize
    dimsg_(1)=in1d%blocksize*in1d%nblocks
    offset_(1)=in1d%offset
    ! if allocated, write the dcontent
    print *, 'dataset=', dataset_id
    if (allocated(in1d%dcontent)) then
      call phdf5_write(in1d%dcontent(1),dims_,dimsg_,offset_,dataset_id)
    elseif (allocated(in1d%zcontent)) then
      call phdf5_write(in1d%zcontent(1),dims_,dimsg_,offset_,dataset_id)
    end if
  end subroutine

!-----------------------------------------------------------------------------
  subroutine get_block1d(in1d,dataset_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block1d), intent(inout) :: in1d
    integer(hid_t) :: dataset_id
    ! local variables
    integer, dimension(1) :: dims_, dimsg_, offset_
    ! set dimension & offset
    dims_(1)=in1d%blocksize
    dimsg_(1)=in1d%blocksize*in1d%nblocks
    offset_(1)=in1d%offset
    ! if allocated, write the dcontent
    if (allocated(in1d%dcontent)) then
      call phdf5_read(in1d%dcontent(1),dims_,dimsg_,offset_,dataset_id)
    elseif (allocated(in1d%zcontent)) then
      call phdf5_read(in1d%zcontent(1),dims_,dimsg_,offset_,dataset_id)
    end if
  end subroutine
end module mod_blocks

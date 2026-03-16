module array_mod
  use, intrinsic :: iso_fortran_env, only : real64
  use MOM_error_infra, only : MOM_err, FATAL
  use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
  implicit none
  private
  public :: RealArray_t, IntArray_t


  type :: RealArray_t
     real(kind=real64), pointer, contiguous :: data(:) => null() !< Storage ptr for array container
     integer :: rank = 0                             !< Rank of array
     integer, allocatable :: shape(:)                !< Shape of array
     integer, allocatable :: lb(:)                   !< Lower bounds
     integer, allocatable :: ub(:)                   !< Upper bounds
   contains
     procedure :: allocReal, freeReal              !< Allocate and deallocate memory in container
     procedure :: viewReal1D,  viewReal2D,  &      !< Associate a Fortran pointer to array container
                  viewReal3D,  viewReal4D
     procedure :: allocReal1D, allocReal2D, &      !< allocate memory and associate a fortran pointer
                  allocReal3D, allocReal4D
     generic   :: view => viewReal1D, viewReal2D, &
                  viewReal3D, viewReal4D           !< Generic interface for view
     generic   :: free => freeReal                 !< Generic interface for deallocate
     generic   :: alloc => allocReal1D, allocReal2D,  &
                  allocReal3D, allocReal4D !< Generic interface for array container allocation
  end type RealArray_t

  type :: IntArray_t
     integer, pointer, contiguous :: data(:) => null() !< Storage ptr for array container
     integer :: rank = 0                   !< Rank of array
     integer, allocatable :: shape(:)      !< Shape of array
     integer, allocatable :: lb(:)         !< Lower bounds
     integer, allocatable :: ub(:)         !< Upper bounds
   contains
     procedure :: allocInt, freeInt                !< Allocate and deallocate memory in container
     procedure ::  viewInt1D,  viewInt2D,  &
                    viewInt3D,  viewInt4D          !< Associate a Fortran pointer to array container
     procedure :: allocInt1D, allocInt2D,  &
                  allocInt3D, allocInt4D          !< Allocate memory and associate a Fortran pointer
     generic   :: view  => viewInt1D, viewInt2D, &
                  viewInt3D, viewInt4D     !< Generic interface for view
     generic   :: free => freeInt          !< Generic interface for deallocate
     generic   :: alloc => allocInt1D, allocInt2D, &
                   allocInt3D, allocInt4D !< Generic interface for array container allocation
  end type intArray_t

contains

subroutine allocReal(this, dims,lb,ub,source)
  class(RealArray_t), intent(inout) :: this         !< The array container to allocate
  integer, intent(in),optional :: dims(:)           !< Dimensions (1-indexed)
  integer, intent(in),optional :: lb(:)             !< Lower bounds
  integer, intent(in),optional :: ub(:)             !< Upper bounds
  real(kind=real64), intent(in), optional :: source !< Initial value for all elements

  if (associated(this%data)) call amrex_deallocate(this%data)
  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)

  if(present(ub) .and. present(lb) .and. .not. present(dims)) then
    if(size(lb) .ne. size(ub)) then
        call MOM_err(FATAL, "allocReal: size of lb and ub must match")
    endif
    this%rank     = size(lb)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = lb(:)
    this%ub(:)    = ub(:)
    this%shape(:) = ub(:)-lb(:)+1
  elseif(present(dims) .and. .not. present(ub) .and. .not. present(lb)) then
    this%rank     = size(dims)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = 1
    this%ub(:)    = dims(:)
    this%shape(:) = dims(:)
  else
    call MOM_err(FATAL, "allocReal: Must specify either ub and lb or dims")
  endif

  ! allocate the memory
  call amrex_allocate(this%data,1,product(this%shape))

  ! initialize the variable
  ! Note this this is a CPU only assignment.
  ! It will not work correctly on the GPU
  if(present(source)) this%data(:) = source


end subroutine allocReal

subroutine allocInt(this, dims,lb,ub,source)
  class(IntArray_t), intent(inout) :: this !< The array container to allocate
  integer, intent(in),optional :: dims(:)  !< Dimensions (1-indexed)
  integer, intent(in),optional :: lb(:)    !< Lower bounds
  integer, intent(in),optional :: ub(:)    !< Upper bounds
  integer, optional :: source              !< Initial value for all elements

  integer :: len                           !< the length of the array to allocate
  integer :: i

  if (associated(this%data)) call amrex_deallocate(this%data)
  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)

  if(present(ub) .and. present(lb) .and. .not. present(dims)) then
    if(size(lb) .ne. size(ub)) then
        call MOM_err(FATAL, "allocReal: size of lb and ub must match")
    endif
    this%rank     = size(lb)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = lb(:)
    this%ub(:)    = ub(:)
    this%shape(:) = ub(:)-lb(:)+1
  elseif(present(dims) .and. .not. present(ub) .and. .not. present(lb)) then
    this%rank     = size(dims)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = 1
    this%ub(:)    = dims(:)
    this%shape(:) = dims(:)
  else
    call MOM_err(FATAL, "allocReal: Must specify either ub and lb or dims")
  endif

  ! allocate the memory
  call amrex_allocate(this%data,1,product(this%shape))

  ! initialize the variable
  ! Note this this is a CPU only assignment.
  ! It will not work correctly on the GPU
  if(present(source)) this%data(:)=source

end subroutine allocInt

subroutine freeReal(this)
  class(RealArray_t), intent(inout) :: this  !< The array container to deallocate

  if (associated(this%data)) call amrex_deallocate(this%data)
  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeReal

subroutine freeInt(this)
  class(IntArray_t), intent(inout) :: this  !< The array container to deallocate

  if (associated(this%data))  call amrex_deallocate(this%data)
  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeInt

subroutine allocReal1D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this         !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:) !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)           !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)             !< Lower bounds
   integer, intent(in),optional :: ub(:)             !< Upper bounds
   real(kind=real64), intent(in), optional :: source !< Initial value for all elements

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine allocReal1D

subroutine viewReal1D(this, a)
   class(RealArray_t), intent(in) :: this     !< The already allocated array container
   real(kind=real64), pointer :: a(:)         !< The Fortran pointer array to associate

   if (this%rank /= 1) call MOM_err(FATAL, "viewReal1D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewReal1D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewReal1D

subroutine allocReal2D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this           !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:) !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)             !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)               !< Lower bounds
   integer, intent(in),optional :: ub(:)               !< Upper bounds
   real(kind=real64), intent(in), optional :: source   !< Initial value for all elements

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine allocReal2D

subroutine viewReal2D(this,a)
   class(RealArray_t), intent(in) :: this              !< The already allocated array container
   real(kind=real64), intent(inout), pointer :: a(:,:) !< The Fortran pointer array to associate

   if (this%rank /= 2) call MOM_err(FATAL, "viewReal2D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewReal2D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine viewReal2D

subroutine allocReal3D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this             !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:,:) !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)               !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                 !< Lower bounds
   integer, intent(in),optional :: ub(:)                 !< Upper bounds
   real(kind=real64), intent(in), optional :: source     !< Initial value for all elements

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine allocReal3D

subroutine viewReal3D(this,a)
   class(RealArray_t), intent(in) :: this                !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:,:) !< The Fortran pointer array

   if (this%rank /= 3) call MOM_err(FATAL, "viewReal3D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewReal3D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine viewReal3D

subroutine allocReal4D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this               !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:,:,:) !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)                 !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                   !< Lower bounds
   integer, intent(in),optional :: ub(:)                   !< Upper bounds
   real(kind=real64), intent(in), optional :: source       !< Initial value for all elements

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine allocReal4D

subroutine viewReal4D(this,a)
   class(RealArray_t), intent(in) :: this                  !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:,:,:) !< The Fortran pointer array

   if (this%rank /= 4) call MOM_err(FATAL, "viewReal4D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewReal4D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine viewReal4D

subroutine allocInt1D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this !< The array container to allocate
   integer, intent(inout), pointer :: a(:)  !< The Fortran pointer array
   integer, intent(in), optional :: dims(:) !< Dimensions (1-indexed)
   integer, intent(in), optional :: lb(:)   !< Lower bounds
   integer, intent(in), optional :: ub(:)   !< Upper bounds
   integer, intent(in), optional :: source  !< Initial value foer all elements

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine allocInt1D

subroutine viewInt1D(this, a)
   class(intArray_t), intent(in) :: this   !< The array container to allocate
   integer, intent(inout), pointer :: a(:) !< The Fortran pointer array

   if (this%rank /= 1) call MOM_err(FATAL, "viewInt1D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewInt1D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewInt1D

subroutine allocInt2D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this  !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:) !< The Fortran pointer array
   integer, intent(in), optional :: dims(:)  !< Dimensions (1-indexed)
   integer, intent(in), optional :: lb(:)    !< Lower bounds
   integer, intent(in), optional :: ub(:)    !< Upper bounds
   integer, intent(in), optional :: source   !< Initial value foer all elements

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine allocInt2D

subroutine viewInt2D(this,a)
   class(intArray_t), intent(in) :: this     !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:) !< The Fortran pointer array

   if (this%rank /= 2) call MOM_err(FATAL, "viewInt2D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewInt2D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine viewInt2D

subroutine allocInt3D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this    !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:,:) !< The Fortran pointer array
   integer, intent(in), optional :: dims(:)    !< Dimensions (1-indexed)
   integer, intent(in), optional :: lb(:)      !< Lower bounds
   integer, intent(in), optional :: ub(:)      !< Upper bounds
   integer, intent(in), optional :: source     !< Initial value foer all elements

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine allocInt3D

subroutine viewInt3D(this,a)
   class(intArray_t), intent(in) :: this       !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:,:) !< The Fortran pointer array

   if (this%rank /= 3) call MOM_err(FATAL, "viewInt3D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewInt3D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine viewInt3D

subroutine allocInt4D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this      !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:,:,:) !< The Fortran pointer array
   integer, intent(in), optional :: dims(:)      !< Dimensions (1-indexed)
   integer, intent(in), optional :: lb(:)        !< Lower bounds
   integer, intent(in), optional :: ub(:)        !< Upper bounds
   integer, intent(in), optional :: source       !< Initial value foer all elements

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine allocInt4D

subroutine viewInt4D(this,a)
   class(intArray_t), intent(in) :: this         !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:,:,:) !< The Fortran pointer array

   if (this%rank /= 4) call MOM_err(FATAL, "viewInt4D: rank mismatch")
   if (.not. allocated(this%shape)) call MOM_err(FATAL, "viewInt4D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine viewInt4D

end module array_mod

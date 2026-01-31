module array_mod
  implicit none
  private
  public :: RealArray_t, IntArray_t, dp

  integer, parameter :: dp = kind(1.0d0)

  type :: RealArray_t
     real(kind=dp), pointer :: data(:)
     integer :: rank = 0
     integer, allocatable :: shape(:)
     integer, allocatable :: lb(:)
     integer, allocatable :: ub(:)
   contains
     procedure :: allocReal, freeReal
     procedure ::  viewReal1D,  viewReal2D,  viewReal3D,  viewReal4D
     procedure :: allocReal1D, allocReal2D, allocReal3D, allocReal4D
     generic   :: view => viewReal1D, viewReal2D, viewReal3D, viewReal4D
     generic   :: free => freeReal
     generic   :: alloc => allocReal1D, allocReal2D, allocReal3D, allocReal4D
     !procedure :: slice_3d
  end type RealArray_t

  type :: IntArray_t
     integer, pointer :: data(:)
     integer :: rank = 0
     integer, allocatable :: shape(:)
     integer, allocatable :: lb(:)
     integer, allocatable :: ub(:)
   contains
     procedure :: allocInt, freeInt
     procedure ::  viewInt1D,  viewInt2D,  viewInt3D,  viewInt4D
     procedure :: allocInt1D, allocInt2D, allocInt3D, allocInt4D
     generic   :: view  => viewInt1D, viewInt2D, viewInt3D, viewInt4D
     generic   :: free => freeInt
     generic   :: alloc => allocInt1D, allocInt2D, allocInt3D, allocInt4D
  end type intArray_t

contains

subroutine allocReal(this, dims,lb,ub,source)
  class(RealArray_t), intent(inout) :: this
  integer, intent(in),optional :: dims(:), lb(:), ub(:)
  real(kind=dp), intent(in), optional :: source

  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)

  if(present(ub) .and. present(lb) .and. .not. present(dims)) then 
    if(size(lb) .ne. size(ub)) then 
        write(*,*) 'Error: size of lb and ub must match'
        !stat=1
        return
    endif
    this%rank     = size(lb)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = lb(:)
    this%ub(:)    = ub(:)
    this%shape(:) = ub(:)-lb(:)+1
    allocate(this%data(product(this%shape)))
    if(present(source)) this%data(:)=source
    !stat=0
  elseif(present(dims) .and. .not. present(ub) .and. .not. present(lb)) then 
    this%rank     = size(dims)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = 1
    this%ub(:)    = dims(:)
    this%shape(:) = dims(:)
    allocate(this%data(product(dims)))
    if(present(source)) this%data(:)=source
    !stat=0
  else
     write(*,*) 'Error: Must specify either ub and lb or dims'
     !stat=1
     return
  endif

end subroutine allocReal

subroutine allocInt(this, dims,lb,ub,source)
  class(IntArray_t), intent(inout) :: this
  integer, intent(in),optional :: dims(:), lb(:), ub(:)
  integer, optional :: source
  !integer, optional :: stat

  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)

  if(present(ub) .and. present(lb) .and. .not. present(dims)) then 
    if(size(lb) .ne. size(ub)) then 
        write(*,*) 'Error: size of lb and ub must match'
        !if(present(stat)) stat=1
        return
    endif
    this%rank     = size(lb)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = lb(:)
    this%ub(:)    = ub(:)
    this%shape(:) = ub(:)-lb(:)+1
    allocate(this%data(product(this%shape)))
    if(present(source)) this%data(:)=source
    !if(present(stat)) stat=0
  elseif(present(dims) .and. .not. present(ub) .and. .not. present(lb)) then 
    this%rank     = size(dims)
    ! Allocate shape and bound information
    allocate(this%shape(this%rank),this%lb(this%rank),this%ub(this%rank))

    this%lb(:)    = 1
    this%ub(:)    = dims(:)
    this%shape(:) = dims(:)
    allocate(this%data(product(dims)))
    if(present(source)) this%data(:)=source
    !if(present(stat)) stat=0
  else
     write(*,*) 'Error: Must specify either ub and lb or dims'
     !if(present(stat)) stat=1
     return
  endif

end subroutine allocInt

subroutine freeReal(this)
  class(RealArray_t), intent(inout) :: this

  if (associated(this%data)) deallocate(this%data)
  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeReal

subroutine freeInt(this)
  class(IntArray_t), intent(inout) :: this

  if (associated(this%data))  deallocate(this%data)
  if (allocated(this%shape)) deallocate(this%shape)
  if (allocated(this%lb))    deallocate(this%lb)
  if (allocated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeInt

subroutine slice_3d(this, sub, i, j, k)
  class(RealArray_t), intent(in) :: this
  real(kind=dp), pointer :: sub(:,:)
  integer, intent(in) :: i(:), j(:), k(:)
  real(kind=dp), pointer :: a(:,:,:)
  integer :: n1, n2, n3
  integer :: iend,jend

  if (this%rank /= 3) stop "slice_3d: rank mismatch"
  ! a => reshape(this%data, this%shape)
  ! Extract dimensions
  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  a(1:n1,1:n2,1:n3) => this%data

  iend=size(i)
  jend=size(j)
  if(is_contiguous(a(i(1):i(iend), j(1):j(jend), k(1)))) then 
     ! Safe for pointr assignment
     sub => a(i(1):i(iend), j(1):j(jend), k(1))
   else
     ! must allocate and copy 
     allocate(sub(iend,jend))
     sub = a(i,j,k(1))
   endif

end subroutine slice_3d

subroutine allocReal1D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this
   real(kind=dp), pointer :: a(:)
   integer, intent(in),optional :: dims(:), lb(:), ub(:)
   real(kind=dp), intent(in), optional :: source

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine allocReal1D

subroutine viewReal1D(this, a)
   class(RealArray_t), intent(in) :: this
   real(kind=dp), pointer :: a(:)
   integer :: n1

   if (this%rank /= 1) stop "viewReal1D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewReal1D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewReal1D

subroutine allocReal2D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this
   real(kind=dp), pointer :: a(:,:)
   integer, intent(in),optional :: dims(:), lb(:), ub(:)
   real(kind=dp), intent(in), optional :: source

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2)) => this%data

end subroutine allocReal2D

subroutine viewReal2D(this,a)
   class(RealArray_t), intent(in) :: this
   real(kind=dp), pointer :: a(:,:)

   if (this%rank /= 2) stop "viewReal2D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewReal2D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2)) => this%data

end subroutine viewReal2D

subroutine allocReal3D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this
   real(kind=dp), pointer :: a(:,:,:)
   integer, intent(in),optional :: dims(:), lb(:), ub(:)
   real(kind=dp), intent(in), optional :: source

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3)) => this%data

end subroutine allocReal3D

subroutine viewReal3D(this,a)
   class(RealArray_t), intent(in) :: this
   real(kind=dp), pointer :: a(:,:,:)

   if (this%rank /= 3) stop "viewReal3D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewReal3D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3)) => this%data

end subroutine viewReal3D

subroutine allocReal4D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this
   real(kind=dp), pointer :: a(:,:,:,:)
   integer, intent(in),optional :: dims(:), lb(:), ub(:)
   real(kind=dp), intent(in), optional :: source

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3),this%lb(4):this%ub(4)) => this%data

end subroutine allocReal4D

subroutine viewReal4D(this,a)
   class(RealArray_t), intent(in) :: this
   real(kind=dp), pointer :: a(:,:,:,:)

   if (this%rank /= 4) stop "viewReal4D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewReal4D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3),this%lb(4):this%ub(4)) => this%data

end subroutine viewReal4D

subroutine allocInt1D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this
   integer, pointer :: a(:)
   integer, intent(in), optional :: dims(:), lb(:), ub(:)
   integer, intent(in), optional :: source

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine allocInt1D

subroutine viewInt1D(this, a)
   class(intArray_t), intent(in) :: this
   integer, pointer :: a(:)

   if (this%rank /= 1) stop "viewInt1D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewInt1D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewInt1D

subroutine allocInt2D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this
   integer, pointer :: a(:,:)
   integer, intent(in), optional :: dims(:), lb(:), ub(:)
   integer, intent(in), optional :: source

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2)) => this%data

end subroutine allocInt2D

subroutine viewInt2D(this,a)
   class(intArray_t), intent(in) :: this
   integer, pointer :: a(:,:)

   if (this%rank /= 2) stop "viewInt2D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewInt2D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2)) => this%data

end subroutine viewInt2D

subroutine allocInt3D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this
   integer, pointer :: a(:,:,:)
   integer, intent(in), optional :: dims(:), lb(:), ub(:)
   integer, intent(in), optional :: source

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3)) => this%data

end subroutine allocInt3D

subroutine viewInt3D(this,a)
   class(intArray_t), intent(in) :: this
   integer, pointer :: a(:,:,:)

   if (this%rank /= 3) stop "viewInt3D: rank mismatch"
   if (.not. allocated(this%shape)) stop "viewInt3D: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3)) => this%data

end subroutine viewInt3D

subroutine allocInt4D(this, a, dims, lb, ub, source)
   class(intArray_t), intent(inout) :: this
   integer, pointer :: a(:,:,:,:)
   integer, intent(in), optional :: dims(:), lb(:), ub(:)
   integer, intent(in), optional :: source

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3),this%lb(4):this%ub(4)) => this%data

end subroutine allocInt4D

subroutine viewInt4D(this,a)
   class(intArray_t), intent(in) :: this
   integer, pointer :: a(:,:,:,:)

   if (this%rank /= 4) stop "view4DInt: rank mismatch"
   if (.not. allocated(this%shape)) stop "view4DInt: shape not allocated"

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1),this%lb(2):this%ub(2),this%lb(3):this%ub(3),this%lb(4):this%ub(4)) => this%data

end subroutine viewInt4D

end module array_mod

module array_mod
  use, intrinsic :: iso_fortran_env, only : real64
  use iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use MOM_error_infra, only : MOM_err, FATAL
  implicit none
  private
  public :: RealArray_t, RealArray_C
  public :: IntArray_t
  public :: LogicalArray_t, LogicalArray_C

  !< Type IntArray_C struct for C++ bridge layer
  type, bind(C) :: IntArray_C
     type(c_ptr) :: data               !< Storage pointer for array container
     type(c_ptr) :: shape              !< An array of dinmension extents
     type(c_ptr) :: lb                 !< Lower bounds
     type(c_ptr) :: ub                 !< Upper bounds
     integer(c_int) :: rank            !< The number of dimensions
  end type IntArray_C

  !< Type RealArray_C struct for C++ bridge layer
  type, bind(C) :: RealArray_C
     type(c_ptr) :: data               !< Storage pointer for array container
     type(c_ptr) :: shape              !< An array of dinmension extents
     type(c_ptr) :: lb                 !< Lower bounds
     type(c_ptr) :: ub                 !< Upper bounds
     integer(c_int) :: rank            !< The number of dimensions
  end type RealArray_C

  !< Type LogicalArray_C struct for C++ bridge layer. The data pointer is
  !! integer-encoded (0/1) -- see LogicalArray_t%to_c in the TIM infra layer.
  type, bind(C) :: LogicalArray_C
     type(c_ptr) :: data               !< Storage pointer for array container
     type(c_ptr) :: shape              !< An array of dinmension extents
     type(c_ptr) :: lb                 !< Lower bounds
     type(c_ptr) :: ub                 !< Upper bounds
     integer(c_int) :: rank            !< The number of dimensions
  end type LogicalArray_C

  type :: RealArray_t
     real(kind=real64), pointer, contiguous :: data(:) => null() !< Storage ptr for array container
     integer :: rank = 0                            !< The number of dimensions
     integer, pointer :: shape(:) => null()         !< An array of dimension extents
     integer, pointer :: lb(:) => null()            !< Lower bounds
     integer, pointer :: ub(:) => null()            !< Upper bounds
   contains
     procedure :: allocReal                    !< Allocate memory in container
     procedure :: freeReal                     !< Deallocates memory from a container
     procedure :: associated => isAssociatedReal !< True if the container holds allocated data
     procedure :: viewReal1D, viewReal2D, &    !< Associates a Fortran pointer to an array container
                  viewReal3D, viewReal4D
     procedure :: allocReal1D, allocReal2D, &  !< Allocate memory and associate a fortran pointer
                  allocReal3D, allocReal4D
     procedure :: copy2FReal1D, copy2FReal2D, & !< Copy data in a RealArray_t to a Fortran array
                  copy2FReal3D, copy2FReal4D
     procedure :: copy2AReal1D, copy2AReal2D, & !< Copy data from a Fortran array to a container
                  copy2AReal3D, copy2AReal4D, &
                  copy2AReal0D
     procedure :: allocViewReal1D, allocViewReal2D, &
                  allocViewReal3D, allocViewReal4D
     procedure :: write_binary                  !< write variable to a binary file
     procedure :: read_binary                   !< read variable from a binary file
     generic :: copy2F => copy2FReal1D, &       !< Generic interface for copy to Fortran arrayc
                copy2FReal2D, copy2FReal3d, &
                copy2FReal4D
     generic :: copy2Array => copy2AReal0D, &   !< Generic interface for copy to array container
                copy2AReal1D, copy2AReal2D, &
                copy2AReal3D, copy2AReal4D
     generic :: view => viewReal1D, &         !< Generic interface for view
                viewReal2D, viewReal3D, &
                viewReal4D
     generic :: alloc => allocReal,       &       !< Generic interface for array container allocation
                allocReal1D, allocReal2D, &
                allocReal3D, allocReal4D
     generic :: allocView =>  &           !< Generic interface for array container allocation and view
                allocViewReal1D, allocViewReal2D, &
                allocViewReal3D, allocViewReal4D
     generic :: free => freeReal              !< Generic interface for deallocate
  end type RealArray_t

  type :: IntArray_t
     integer, pointer, contiguous :: data(:) => null() !< Storage ptr for array container
     integer :: rank = 0                     !< The number of dimensions
     integer, pointer :: shape(:) => null()  !< An array of dimension extents
     integer, pointer :: lb(:) => null()     !< Lower bounds
     integer, pointer :: ub(:) => null()     !< Upper bounds
   contains
     procedure :: allocInt                   !< Allocates  memory in container
     procedure :: freeInt                    !< Deallocates memory from a container
     procedure :: associated => isAssociatedInt !< True if the container holds allocated data
     procedure ::  viewInt1D,  viewInt2D, &   !< Associates a Fortran pointer to an array container
                   viewInt3D,  viewInt4D
    procedure :: allocInt1D, allocInt2D,  &  !< Allocates memory and optionally initializes
                  allocInt3D, allocInt4D
     procedure :: allocViewInt1D, allocViewInt2D, & !< Allocates memory and associates a Fortran pointer
                  allocViewInt3D, allocViewInt4D
     procedure :: copy2FInt1D, copy2FInt2D, & !< Copy data in a IntArray_t to a Fortran array
                  copy2FInt3D, copy2FInt4D
     procedure :: copy2AInt1D, copy2AInt2D, & !< Copy data from a Fortran array to an  container
                  copy2AInt3D, copy2AInt4D, &
                  copy2AInt0D
     generic :: copy2F =>                 &       !< Generic interface for copy to Fortran arrayc
                copy2FInt1D, copy2FInt2D, &
                copy2FInt3D, copy2Fint4D
     generic :: copy2Array => copy2AInt0D, &   !< Generic interface for copy to array container
                copy2AInt1D, copy2AInt2D, &
                copy2AInt3D, copy2AInt4D
     generic   :: view  => viewInt1D, &             !< Generic interface for view
                  viewInt2D, viewInt3D, viewInt4D
     generic   :: alloc => allocInt, &      !< Generic interface for array container allocation
                  allocInt1D, allocInt2D, &
                  allocInt3D, allocInt4D
     generic :: allocView =>   &   !< Generic interface for array container allocation and view
                allocViewInt1D, allocViewInt2D, &
                allocViewInt3D, allocViewInt4D
     generic   :: free => freeInt          !< Generic interface for deallocate
  end type intArray_t

  type :: LogicalArray_t
     logical, pointer, contiguous :: data(:) => null() !< Storage ptr for array container
     integer :: rank = 0                     !< The number of dimensions
     integer, pointer :: shape(:) => null()  !< An array of dimension extents
     integer, pointer :: lb(:) => null()     !< Lower bounds
     integer, pointer :: ub(:) => null()     !< Upper bounds
   contains
     procedure :: allocLogical                       !< Allocates memory in container
     procedure :: freeLogical                        !< Deallocates memory from a container
     procedure :: associated => isAssociatedLogical  !< True if the container holds allocated data
     procedure ::  viewLogical1D,  viewLogical2D, &   !< Associates a Fortran pointer to an array container
                   viewLogical3D,  viewLogical4D
     procedure :: allocLogical1D, allocLogical2D,  &  !< Allocates memory and associates a Fortran pointer
                  allocLogical3D, allocLogical4D
     procedure :: allocViewLogical1D, allocViewLogical2D, & !< Allocates memory and associates a Fortran pointer
                  allocViewLogical3D, allocViewLogical4D
     procedure :: copy2FLogical1D, copy2FLogical2D, & !< Copy data in a LogicalArray_t to a Fortran array
                  copy2FLogical3D, copy2FLogical4D
     procedure :: copy2ALogical1D, copy2ALogical2D, & !< Copy data from a Fortran array to a container
                  copy2ALogical3D, copy2ALogical4D, &
                  copy2ALogical0D
     generic :: copy2F =>                 &       !< Generic interface for copy to Fortran arrayc
                copy2FLogical1D, copy2FLogical2D, &
                copy2FLogical3D, copy2FLogical4D
     generic :: copy2Array => copy2ALogical0D, &   !< Generic interface for copy to array container
                copy2ALogical1D, copy2ALogical2D, &
                copy2ALogical3D, copy2ALogical4D
     generic   :: view  => viewLogical1D, &             !< Generic interface for view
                  viewLogical2D, viewLogical3D, viewLogical4D
     generic   :: alloc => allocLogical, &      !< Generic interface for array container allocation
                  allocLogical1D, allocLogical2D, &
                  allocLogical3D, allocLogical4D
     generic :: allocView =>   &   !< Generic interface for array container allocation and view
                allocViewLogical1D, allocViewLogical2D, &
                allocViewLogical3D, allocViewLogical4D
     generic   :: free => freeLogical          !< Generic interface for deallocate
  end type LogicalArray_t

contains

!< This subroutine writes an RealArray_t variable to a binary file
subroutine write_binary(this, unit)
  class(RealArray_t), intent(in) :: this  !< The RealArray_t variable to write to disk
  integer,            intent(in) :: unit  !< The file unit

  ! local variables
  integer :: i
  integer :: n

  ! --- Null case ---
  if (.not. associated(this%data)) then
    write(unit) -1   ! rank = -1 signals null
    return
  endif

  ! --- Rank ---
  n = this%rank
  write(unit) n

  ! --- Write shape ---
  do i=1,n
    write(unit) this%shape(i)
  enddo

  ! --- Write bounds ---
  do i=1,n
    write(unit) this%lb(i)
    write(unit) this%ub(i)
  enddo

  ! --- Write data size ---
  write(unit) size(this%data)

  ! --- Write payload ---
  write(unit) this%data

end subroutine write_binary

!< This subroutine reads an RealArray_t variable from a binary file
subroutine read_binary(this, unit)
  class(RealArray_t), intent(inout) :: this  !< The RealArray_t variable to read from disk
  integer,            intent(in)    :: unit  !< the file unit

  ! local variables
  integer :: i
  integer :: n
  integer :: total_size

  ! --- Read rank ---
  read(unit) n

  ! --- Null case ---
  if (n == -1) then
    if (associated(this%data)) deallocate(this%data)
    if (associated(this%shape)) deallocate(this%shape)
    if (associated(this%lb)) deallocate(this%lb)
    if (associated(this%ub)) deallocate(this%ub)

    nullify(this%data)
    nullify(this%shape)
    nullify(this%lb)
    nullify(this%ub)
    this%rank = 0
    return
  endif

  this%rank = n

  ! --- Clean old allocations ---
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb)) deallocate(this%lb)
  if (associated(this%ub)) deallocate(this%ub)
  if (associated(this%data)) deallocate(this%data)

  ! --- Allocate metadata ---
  allocate(this%shape(n))
  allocate(this%lb(n))
  allocate(this%ub(n))

  ! --- Read shape ---
  do i=1,n
    read(unit) this%shape(i)
  enddo

  ! --- Read bounds ---
  do i=1,n
    read(unit) this%lb(i)
    read(unit) this%ub(i)
  enddo

  ! --- Read data size ---
  read(unit) total_size

  ! --- Allocate and read data ---
  if (total_size > 0) then
    allocate(this%data(total_size))
    read(unit) this%data
  else
    nullify(this%data)
  endif

end subroutine read_binary

!< Function to convert a Fortran structure to a C structure
!function to_c_Real(this) result(cdesc)
!  class(RealArray_t), intent(in) :: this  !< RealArray_t structure to convert to C
!  type(RealArray_C) :: cdesc              !< Resulting C structure
!
!  cdesc%data  = c_loc(this%data(1))
!  cdesc%shape = c_loc(this%shape(1))
!  cdesc%lb    = c_loc(this%lb(1))
!  cdesc%ub    = c_loc(this%ub(1))
!  cdesc%rank  = this%rank
!end function to_c_Real
!
!!< Function to convert a Fortran structure to a C structure
!function to_c_Int(this) result(cdesc)
!  class(IntArray_t), intent(in) :: this    !< IntArray_t structure to convert to C
!  type(IntArray_C) :: cdesc                !< Resulting C structure
!
!  cdesc%data  = c_loc(this%data(1))
!  cdesc%shape = c_loc(this%shape(1))
!  cdesc%lb    = c_loc(this%lb(1))
!  cdesc%ub    = c_loc(this%ub(1))
!  cdesc%rank  = this%rank
!end function to_c_Int
!
!!< Function to convert a Fortran structure to a C structure. this%data is
!!! logical; cdesc%data must point at an integer-encoded (0/1) shadow
!!! buffer, not this%data itself -- see the TIM infra layer for the
!!! working implementation.
!function to_c_Logical(this) result(cdesc)
!  class(LogicalArray_t), intent(inout) :: this !< LogicalArray_t structure to convert to C
!  type(LogicalArray_C) :: cdesc                !< Resulting C structure
!
!  cdesc%data  = c_loc(this%data_c(1))
!  cdesc%shape = c_loc(this%shape(1))
!  cdesc%lb    = c_loc(this%lb(1))
!  cdesc%ub    = c_loc(this%ub(1))
!  cdesc%rank  = this%rank
!end function to_c_Logical

subroutine allocReal(this, dims,lb,ub,source)
  class(RealArray_t), intent(inout) :: this         !< The array container to allocate
  integer, intent(in),optional :: dims(:)           !< Dimensions (1-indexed)
  integer, intent(in),optional :: lb(:)             !< Lower bounds
  integer, intent(in),optional :: ub(:)             !< Upper bounds
  real(kind=real64), intent(in), optional :: source !< Initial value for all elements

  if (associated(this%data))  deallocate(this%data)
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb))    deallocate(this%lb)
  if (associated(this%ub))    deallocate(this%ub)

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
  allocate(this%data(product(this%shape)))

  ! initialize the variable
  ! Note this this is a CPU only assignment.
  ! It will not work correctly on the GPU
  if(present(source)) this%data(:) = source

end subroutine allocReal

subroutine copy2AReal0D(this,var)
  class(RealArray_t), intent(inout) :: this  !< The destination array container
  real, intent(in) :: var  !< The source Fortran array

  ! Local variables
  integer :: i, n

  n  = product(this%shape)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n)
     this%data(i) = var
  enddo

end subroutine copy2AReal0D

!< Copy from 1D Fortran array to RealArray_t
subroutine copy2AReal1D(this,var)
  class(RealArray_t), intent(inout) :: this  !< The destination array container
  real, dimension(:), intent(in) :: var      !< The source Fortran array

  ! Local variables
  integer :: i, n1

  n1 = this%shape(1)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     this%data(i) = var(i)
  enddo

end subroutine copy2AReal1D

!< Copy from 2D Fortran array to RealArray_t
subroutine copy2AReal2D(this,var)
  class(RealArray_t), intent(inout) :: this  !< The destination array container
  real, dimension(:,:), intent(in) :: var    !< The source Fortran array

  ! Local variables
  integer :: i, j, n1, n2

  n1 = this%shape(1)
  n2 = this%shape(2)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (j=1:n2,i=1:n1)
     this%data(i + n1*(j-1)) = var(i,j)
  enddo

end subroutine copy2AReal2D

!< Copy from 3D Fortran array to RealArray_t
subroutine copy2AReal3D(this,var)
  class(RealArray_t), intent(inout) :: this  !< The destination array container
  real, dimension(:,:,:), intent(in) :: var  !< The source Fortran array

  ! Local variables
  integer :: i, j, k, n1, n2, n3

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (k=1:n3, j=1:n2, i=1:n1)
     this%data(i + n1*(j-1) + n1*n2*(k-1)) = var(i,j,k)
  enddo

end subroutine copy2AReal3D

!< Copy from 4D Fortran array to RealArray_t
subroutine copy2AReal4D(this,var)
  class(RealArray_t), intent(inout) :: this   !< The destination array container
  real, dimension(:,:,:,:), intent(in) :: var !< The source Fortran array

  ! Local variables
  integer :: i, j, k, m, n1, n2, n3, n4

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  n4 = this%shape(4)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (m=1:n4, k=1:n3, j=1:n2, i=1:n1)
     this%data(i + n1*(j-1) + n1*n2*(k-1) + n1*n2*n3*(m-1)) = var(i,j,k,m)
  enddo

end subroutine copy2AReal4D

! Copy from 1D RealArray_t to Fortran
subroutine copy2FReal1D(this,var)
  class(RealArray_t), intent(in) :: this    !< The source array container
  real, dimension(:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, n1

  n1 = this%shape(1)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     var(i) = this%data(i)
  enddo

end subroutine copy2FReal1D

! Copy from 2D RealArray_t to Fortran
subroutine copy2FReal2D(this,var)
  class(RealArray_t), intent(in) :: this     !< The source array container
  real, dimension(:,:), intent(inout) :: var !< The destination Fortran array

  ! Local variables
  integer :: i, j, n1, n2

  n1 = this%shape(1)
  n2 = this%shape(2)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (j=1:n2, i=1:n1)
     var(i,j) = this%data(i + n1*(j-1))
  enddo

end subroutine copy2FReal2D

! Copy from 3D RealArray_t to Fortran
subroutine copy2FReal3D(this,var)
  class(RealArray_t), intent(in) :: this       !< The source array container
  real, dimension(:,:,:), intent(inout) :: var !< The destination Fortran array

  ! Local variables
  integer :: i, j, k, n1, n2, n3

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (k=1:n3, j=1:n2, i=1:n1)
     var(i,j,k) = this%data(i + n1*(j-1) + n1*n2*(k-1))
  enddo

end subroutine copy2FReal3D

! Copy from 4D RealArray_t to Fortran
subroutine copy2FReal4D(this,var)
  class(RealArray_t), intent(in) :: this          !< The source array container
  real, dimension(:,:,:,:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, j, k, m, n1, n2, n3, n4

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  n4 = this%shape(4)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (m=1:n4, k=1:n3, j=1:n2, i=1:n1)
     var(i,j,k,m) = this%data(i + n1*(j-1) + n1*n2*(k-1) + n1*n2*n3*(m-1))
  enddo

end subroutine copy2FReal4D

!< Copy from Fortran array to IntArray_t
subroutine copy2AInt0D(this,var)
  class(IntArray_t), intent(inout) :: this  !< The destination array container
  integer, intent(in) :: var      !< The source Fortran array

  ! Local variables
  integer :: i, n1

  n1 = product(this%shape)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     this%data(i) = var
  enddo

end subroutine copy2AInt0D

!< Copy from Fortran array to IntArray_t
subroutine copy2AInt1D(this,var)
  class(IntArray_t), intent(inout) :: this  !< The destination array container
  integer, dimension(:), intent(in) :: var      !< The source Fortran array

  ! Local variables
  integer :: i, n1

  n1 = this%shape(1)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     this%data(i) = var(i)
  enddo

end subroutine copy2AInt1D

!< Copy from Fortran array to IntArray_t
subroutine copy2AInt2D(this,var)
  class(IntArray_t), intent(inout) :: this  !< The destination array container
  integer, dimension(:,:), intent(in) :: var    !< The source Fortran array

  ! Local variables
  integer :: i, j, n1, n2

  n1 = this%shape(1)
  n2 = this%shape(2)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (j=1:n2,i=1:n1)
     this%data(i + n1*(j-1)) = var(i,j)
  enddo

end subroutine copy2AInt2D

!< Copy from Fortran array to IntArray_t
subroutine copy2AInt3D(this,var)
  class(IntArray_t), intent(inout) :: this  !< The destination array container
  integer, dimension(:,:,:), intent(in) :: var  !< The source Fortran array

  ! Local variables
  integer :: i, j, k, n1, n2, n3

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (k=1:n3, j=1:n2, i=1:n1)
     this%data(i + n1*(j-1) + n1*n2*(k-1)) = var(i,j,k)
  enddo

end subroutine copy2AInt3D

!< Copy from Fortran array to IntArray_t
subroutine copy2AInt4D(this,var)
  class(IntArray_t), intent(inout) :: this   !< The destination array container
  integer, dimension(:,:,:,:), intent(in) :: var !< The source Fortran array

  ! Local variables
  integer :: i, j, k, m, n1, n2, n3, n4

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  n4 = this%shape(4)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (m=1:n4, k=1:n3, j=1:n2, i=1:n1)
     this%data(i + n1*(j-1) + n1*n2*(k-1) + n1*n2*n3*(m-1)) = var(i,j,k,m)
  enddo

end subroutine copy2AInt4D

! Copy from 1D IntArray_t to Fortran
subroutine copy2FInt1D(this,var)
  class(IntArray_t), intent(in) :: this    !< The source array container
  integer, dimension(:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, n1

  n1 = this%shape(1)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     var(i) = this%data(i)
  enddo

end subroutine copy2FInt1D

! Copy from 2D IntArray_t to Fortran
subroutine copy2FInt2D(this,var)
  class(IntArray_t), intent(in) :: this     !< The source array container
  integer, dimension(:,:), intent(inout) :: var !< The destination Fortran array

  integer :: i,j, n1, n2
  n1 = this%shape(1)
  n2 = this%shape(2)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (j=1:n2, i=1:n1)
     var(i,j) = this%data(i + n1*(j-1))
  enddo

end subroutine copy2FInt2D

! Copy from 3D IntArray_t to Fortran
subroutine copy2FInt3D(this,var)
  class(IntArray_t), intent(in) :: this       !< The source array container
  integer, dimension(:,:,:), intent(inout) :: var !< The destination Fortran array

  ! Local variables
  integer :: i, j, k, n1, n2, n3

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (k=1:n3, j=1:n2, i=1:n1)
     var(i,j,k) = this%data(i + n1*(j-1) + n1*n2*(k-1))
  enddo

end subroutine copy2FInt3D

! Copy from 4D IntArray_t to Fortran
subroutine copy2FInt4D(this,var)
  class(IntArray_t), intent(in) :: this          !< The source array container
  integer, dimension(:,:,:,:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, j, k, m, n1, n2, n3, n4

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  n4 = this%shape(4)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (m=1:n4, k=1:n3, j=1:n2, i=1:n1)
     var(i,j,k,m) = this%data(i + n1*(j-1) + n1*n2*(k-1) + n1*n2*n3*(m-1))
  enddo

end subroutine copy2FInt4D


subroutine allocInt(this, dims,lb,ub,source)
  class(IntArray_t), intent(inout) :: this !< The array container to allocate
  integer, intent(in),optional :: dims(:)  !< Dimensions (1-indexed)
  integer, intent(in),optional :: lb(:)    !< Lower bounds
  integer, intent(in),optional :: ub(:)    !< Upper bounds
  integer, optional :: source              !< Initial value for all elements

  integer :: len                           !< the length of the array to allocate

  if (associated(this%data)) deallocate(this%data)
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb))    deallocate(this%lb)
  if (associated(this%ub))    deallocate(this%ub)

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
  allocate(this%data(product(this%shape)))

  ! initialize the variable
  ! Note this this is a CPU only assignment.
  ! It will not work correctly on the GPU
  if(present(source)) call this%copy2Array(source)

end subroutine allocInt

subroutine freeReal(this)
  class(RealArray_t), intent(inout) :: this  !< The array container to deallocate

  if (associated(this%data))  deallocate(this%data)
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb))    deallocate(this%lb)
  if (associated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeReal

subroutine freeInt(this)
  class(IntArray_t), intent(inout) :: this  !< The array container to deallocate

  if (associated(this%data))  deallocate(this%data)
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb))    deallocate(this%lb)
  if (associated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeInt

pure function isAssociatedReal(this) result(is_assoc)
  class(RealArray_t), intent(in) :: this  !< The array container to query
  logical :: is_assoc                     !< True if the container holds allocated data

  is_assoc = associated(this%data)
end function isAssociatedReal

pure function isAssociatedInt(this) result(is_assoc)
  class(IntArray_t), intent(in) :: this  !< The array container to query
  logical :: is_assoc                    !< True if the container holds allocated data

  is_assoc = associated(this%data)
end function isAssociatedInt

subroutine allocReal1D(this, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   real(kind=real64), intent(in) :: source(:) !< Assignment array

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub)

   ! copy the array in
   call this%copy2AReal1D(source)

end subroutine allocReal1D

subroutine allocViewReal1D(this, a, dims, lb, ub, source)
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

end subroutine allocViewReal1D

subroutine viewReal1D(this, a)
   class(RealArray_t), intent(in) :: this     !< The already allocated array container
   real(kind=real64), pointer :: a(:)         !< The Fortran pointer array to associate

   if (this%rank /= 1) call MOM_err(FATAL, "viewReal1D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewReal1D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewReal1D

subroutine allocReal2D(this, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   real(kind=real64), intent(in) :: source(:,:) !< Assignment array

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub)

   call this%copy2AReal2D(source)

end subroutine allocReal2D

subroutine allocViewReal2D(this, a, dims, lb, ub, source)
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

end subroutine allocViewReal2D

subroutine viewReal2D(this,a)
   class(RealArray_t), intent(in) :: this              !< The already allocated array container
   real(kind=real64), intent(inout), pointer :: a(:,:) !< The Fortran pointer array to associate

   if (this%rank /= 2) call MOM_err(FATAL, "viewReal2D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewReal2D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine viewReal2D

subroutine allocReal3D(this, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   real(kind=real64), intent(in) :: source(:,:,:) !< Assignment array

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub)

   call this%copy2AReal3D(source)

end subroutine allocReal3D

subroutine allocViewReal3D(this, a, dims, lb, ub, source)
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

end subroutine allocViewReal3D

subroutine viewReal3D(this,a)
   class(RealArray_t), intent(in) :: this                !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:,:) !< The Fortran pointer array

   if (this%rank /= 3) call MOM_err(FATAL, "viewReal3D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewReal3D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine viewReal3D

subroutine allocReal4D(this, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   real(kind=real64), intent(in) :: source(:,:,:,:) !< Assignment array

   ! allocate the memory
   call this%allocReal(dims=dims, lb=lb, ub=ub)

   call this%copy2AReal4D(source)

end subroutine allocReal4D

subroutine allocViewReal4D(this, a, dims, lb, ub, source)
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

end subroutine allocViewReal4D

subroutine viewReal4D(this,a)
   class(RealArray_t), intent(in) :: this                  !< The array container to allocate
   real(kind=real64), intent(inout), pointer :: a(:,:,:,:) !< The Fortran pointer array

   if (this%rank /= 4) call MOM_err(FATAL, "viewReal4D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewReal4D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine viewReal4D

subroutine allocInt1D(this, dims, lb, ub, source)
   class(IntArray_t), intent(inout) :: this    !< The array container to allocate
   integer, intent(in),optional :: dims(:)     !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)       !< Lower bounds
   integer, intent(in),optional :: ub(:)       !< Upper bounds
   integer, intent(in)          :: source(:)   !< Assignment array

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub)

   call this%copy2AInt1D(source)

end subroutine allocInt1D

subroutine allocViewInt1D(this, a, dims, lb, ub, source)
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

end subroutine allocViewInt1D

subroutine viewInt1D(this, a)
   class(intArray_t), intent(in) :: this   !< The array container to allocate
   integer, intent(inout), pointer :: a(:) !< The Fortran pointer array

   if (this%rank /= 1) call MOM_err(FATAL, "viewInt1D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewInt1D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewInt1D

subroutine allocInt2D(this, dims, lb, ub, source)
   class(IntArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   integer, intent(in)          :: source(:,:) !< Assignment array

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub)

   call this%copy2Array(source)

end subroutine allocInt2D

subroutine allocViewInt2D(this, a, dims, lb, ub, source)
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

end subroutine allocViewInt2D

subroutine viewInt2D(this,a)
   class(intArray_t), intent(in) :: this     !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:) !< The Fortran pointer array

   if (this%rank /= 2) call MOM_err(FATAL, "viewInt2D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewInt2D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine viewInt2D

subroutine allocInt3D(this, dims, lb, ub, source)
   class(IntArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   integer, intent(in) :: source(:,:,:) !< Assignment array

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub)

   call this%copy2Array(source)

end subroutine allocInt3D

subroutine allocViewInt3D(this, a, dims, lb, ub, source)
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

end subroutine allocViewInt3D

subroutine viewInt3D(this,a)
   class(intArray_t), intent(in) :: this       !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:,:) !< The Fortran pointer array

   if (this%rank /= 3) call MOM_err(FATAL, "viewInt3D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewInt3D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine viewInt3D

subroutine allocInt4D(this, dims, lb, ub, source)
   class(IntArray_t), intent(inout) :: this                !< The array container to allocate
   integer, intent(in),optional :: dims(:)                  !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)                    !< Lower bounds
   integer, intent(in),optional :: ub(:)                    !< Upper bounds
   integer, intent(in) :: source(:,:,:,:) !< Assignment array

   ! allocate the memory
   call this%allocInt(dims=dims, lb=lb, ub=ub)

   ! assign the values in the array container
   call this%copy2Array(source)

end subroutine allocInt4D

subroutine allocViewInt4D(this, a, dims, lb, ub, source)
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

end subroutine allocViewInt4D

subroutine viewInt4D(this,a)
   class(intArray_t), intent(in) :: this         !< The array container to allocate
   integer, intent(inout), pointer :: a(:,:,:,:) !< The Fortran pointer array

   if (this%rank /= 4) call MOM_err(FATAL, "viewInt4D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewInt4D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine viewInt4D

subroutine allocLogical(this, dims,lb,ub,source)
  class(LogicalArray_t), intent(inout) :: this !< The array container to allocate
  integer, intent(in),optional :: dims(:)      !< Dimensions (1-indexed)
  integer, intent(in),optional :: lb(:)        !< Lower bounds
  integer, intent(in),optional :: ub(:)        !< Upper bounds
  logical, intent(in), optional :: source      !< Initial value for all elements

  if (associated(this%data))  deallocate(this%data)
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb))    deallocate(this%lb)
  if (associated(this%ub))    deallocate(this%ub)

  if(present(ub) .and. present(lb) .and. .not. present(dims)) then
    if(size(lb) .ne. size(ub)) then
        call MOM_err(FATAL, "allocLogical: size of lb and ub must match")
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
    call MOM_err(FATAL, "allocLogical: Must specify either ub and lb or dims")
  endif

  ! allocate the memory
  allocate(this%data(product(this%shape)))

  ! initialize the variable
  ! Note this this is a CPU only assignment.
  ! It will not work correctly on the GPU
  if(present(source)) call this%copy2Array(source)

end subroutine allocLogical

subroutine freeLogical(this)
  class(LogicalArray_t), intent(inout) :: this  !< The array container to deallocate

  if (associated(this%data))  deallocate(this%data)
  if (associated(this%shape)) deallocate(this%shape)
  if (associated(this%lb))    deallocate(this%lb)
  if (associated(this%ub))    deallocate(this%ub)
  this%rank = 0
end subroutine freeLogical

pure function isAssociatedLogical(this) result(is_assoc)
  class(LogicalArray_t), intent(in) :: this  !< The array container to query
  logical :: is_assoc                        !< True if the container holds allocated data

  is_assoc = associated(this%data)
end function isAssociatedLogical

subroutine copy2ALogical0D(this,var)
  class(LogicalArray_t), intent(inout) :: this  !< The destination array container
  logical, intent(in) :: var  !< The source Fortran scalar

  ! Local variables
  integer :: i, n

  n  = product(this%shape)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n)
     this%data(i) = var
  enddo

end subroutine copy2ALogical0D

!< Copy from 1D Fortran array to LogicalArray_t
subroutine copy2ALogical1D(this,var)
  class(LogicalArray_t), intent(inout) :: this  !< The destination array container
  logical, dimension(:), intent(in) :: var      !< The source Fortran array

  ! Local variables
  integer :: i, n1

  n1 = this%shape(1)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     this%data(i) = var(i)
  enddo

end subroutine copy2ALogical1D

!< Copy from 2D Fortran array to LogicalArray_t
subroutine copy2ALogical2D(this,var)
  class(LogicalArray_t), intent(inout) :: this  !< The destination array container
  logical, dimension(:,:), intent(in) :: var    !< The source Fortran array

  ! Local variables
  integer :: i, j, n1, n2

  n1 = this%shape(1)
  n2 = this%shape(2)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (j=1:n2,i=1:n1)
     this%data(i + n1*(j-1)) = var(i,j)
  enddo

end subroutine copy2ALogical2D

!< Copy from 3D Fortran array to LogicalArray_t
subroutine copy2ALogical3D(this,var)
  class(LogicalArray_t), intent(inout) :: this  !< The destination array container
  logical, dimension(:,:,:), intent(in) :: var  !< The source Fortran array

  ! Local variables
  integer :: i, j, k, n1, n2, n3

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (k=1:n3, j=1:n2, i=1:n1)
     this%data(i + n1*(j-1) + n1*n2*(k-1)) = var(i,j,k)
  enddo

end subroutine copy2ALogical3D

!< Copy from 4D Fortran array to LogicalArray_t
subroutine copy2ALogical4D(this,var)
  class(LogicalArray_t), intent(inout) :: this    !< The destination array container
  logical, dimension(:,:,:,:), intent(in) :: var  !< The source Fortran array

  ! Local variables
  integer :: i, j, k, m, n1, n2, n3, n4

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  n4 = this%shape(4)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (m=1:n4, k=1:n3, j=1:n2, i=1:n1)
     this%data(i + n1*(j-1) + n1*n2*(k-1) + n1*n2*n3*(m-1)) = var(i,j,k,m)
  enddo

end subroutine copy2ALogical4D

! Copy from 1D LogicalArray_t to Fortran
subroutine copy2FLogical1D(this,var)
  class(LogicalArray_t), intent(in) :: this    !< The source array container
  logical, dimension(:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, n1

  n1 = this%shape(1)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (i=1:n1)
     var(i) = this%data(i)
  enddo

end subroutine copy2FLogical1D

! Copy from 2D LogicalArray_t to Fortran
subroutine copy2FLogical2D(this,var)
  class(LogicalArray_t), intent(in) :: this      !< The source array container
  logical, dimension(:,:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, j, n1, n2

  n1 = this%shape(1)
  n2 = this%shape(2)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (j=1:n2, i=1:n1)
     var(i,j) = this%data(i + n1*(j-1))
  enddo

end subroutine copy2FLogical2D

! Copy from 3D LogicalArray_t to Fortran
subroutine copy2FLogical3D(this,var)
  class(LogicalArray_t), intent(in) :: this        !< The source array container
  logical, dimension(:,:,:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, j, k, n1, n2, n3

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (k=1:n3, j=1:n2, i=1:n1)
     var(i,j,k) = this%data(i + n1*(j-1) + n1*n2*(k-1))
  enddo

end subroutine copy2FLogical3D

! Copy from 4D LogicalArray_t to Fortran
subroutine copy2FLogical4D(this,var)
  class(LogicalArray_t), intent(in) :: this          !< The source array container
  logical, dimension(:,:,:,:), intent(inout) :: var  !< The destination Fortran array

  ! Local variables
  integer :: i, j, k, m, n1, n2, n3, n4

  n1 = this%shape(1)
  n2 = this%shape(2)
  n3 = this%shape(3)
  n4 = this%shape(4)

  ! do concurrent so the copy runs on the device under offload
  do concurrent (m=1:n4, k=1:n3, j=1:n2, i=1:n1)
     var(i,j,k,m) = this%data(i + n1*(j-1) + n1*n2*(k-1) + n1*n2*n3*(m-1))
  enddo

end subroutine copy2FLogical4D

subroutine allocLogical1D(this, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in) :: source(:)              !< Assignment array

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub)

   ! copy the array in
   call this%copy2ALogical1D(source)

end subroutine allocLogical1D

subroutine allocViewLogical1D(this, a, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   logical, intent(inout), pointer :: a(:)       !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in), optional :: source       !< Initial value for all elements

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine allocViewLogical1D

subroutine viewLogical1D(this, a)
   class(LogicalArray_t), intent(in) :: this  !< The already allocated array container
   logical, pointer :: a(:)                   !< The Fortran pointer array to associate

   if (this%rank /= 1) call MOM_err(FATAL, "viewLogical1D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewLogical1D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1)) => this%data

end subroutine viewLogical1D

subroutine allocLogical2D(this, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in) :: source(:,:)            !< Assignment array

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub)

   call this%copy2ALogical2D(source)

end subroutine allocLogical2D

subroutine allocViewLogical2D(this, a, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   logical, intent(inout), pointer :: a(:,:)     !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in), optional :: source       !< Initial value for all elements

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine allocViewLogical2D

subroutine viewLogical2D(this,a)
   class(LogicalArray_t), intent(in) :: this   !< The already allocated array container
   logical, intent(inout), pointer :: a(:,:)   !< The Fortran pointer array to associate

   if (this%rank /= 2) call MOM_err(FATAL, "viewLogical2D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewLogical2D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2)) => this%data

end subroutine viewLogical2D

subroutine allocLogical3D(this, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in) :: source(:,:,:)          !< Assignment array

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub)

   call this%copy2ALogical3D(source)

end subroutine allocLogical3D

subroutine allocViewLogical3D(this, a, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   logical, intent(inout), pointer :: a(:,:,:)   !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in), optional :: source       !< Initial value for all elements

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine allocViewLogical3D

subroutine viewLogical3D(this,a)
   class(LogicalArray_t), intent(in) :: this    !< The array container to allocate
   logical, intent(inout), pointer :: a(:,:,:)  !< The Fortran pointer array

   if (this%rank /= 3) call MOM_err(FATAL, "viewLogical3D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewLogical3D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3)) => this%data

end subroutine viewLogical3D

subroutine allocLogical4D(this, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this  !< The array container to allocate
   integer, intent(in),optional :: dims(:)       !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)         !< Lower bounds
   integer, intent(in),optional :: ub(:)         !< Upper bounds
   logical, intent(in) :: source(:,:,:,:)        !< Assignment array

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub)

   ! assign the values in the array container
   call this%copy2ALogical4D(source)

end subroutine allocLogical4D

subroutine allocViewLogical4D(this, a, dims, lb, ub, source)
   class(LogicalArray_t), intent(inout) :: this   !< The array container to allocate
   logical, intent(inout), pointer :: a(:,:,:,:)  !< The Fortran pointer array
   integer, intent(in),optional :: dims(:)        !< Dimensions (1-indexed)
   integer, intent(in),optional :: lb(:)          !< Lower bounds
   integer, intent(in),optional :: ub(:)          !< Upper bounds
   logical, intent(in), optional :: source        !< Initial value for all elements

   ! allocate the memory
   call this%allocLogical(dims=dims, lb=lb, ub=ub, source=source)

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine allocViewLogical4D

subroutine viewLogical4D(this,a)
   class(LogicalArray_t), intent(in) :: this      !< The array container to allocate
   logical, intent(inout), pointer :: a(:,:,:,:)  !< The Fortran pointer array

   if (this%rank /= 4) call MOM_err(FATAL, "viewLogical4D: rank mismatch")
   if (.not. associated(this%shape)) call MOM_err(FATAL, "viewLogical4D: shape not allocated")

   ! Zero copy no allocation
   a(this%lb(1):this%ub(1), this%lb(2):this%ub(2), &
     this%lb(3):this%ub(3), this%lb(4):this%ub(4)) => this%data

end subroutine viewLogical4D



end module array_mod

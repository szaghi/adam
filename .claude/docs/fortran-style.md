# Fortran Style Guide

General Fortran coding guidelines for this codebase, inspired by the [Zen of Fortran](https://github.com/szaghi/zen-of-fortran).

## The Zen of Fortran

1. **Standard compliance** is better than *custom or extended*
2. **Beautiful** is better than *ugly*:
   - **Readability** counts
   - **Explicit** is better than *impl.*
   - **Simple** is better than *CoMpleX*
   - **CoMpleX** is better than *c0mp1|c@ted*
   - **Flat** is better than *nested*
   - **S p a r s e** is better than *dense*
3. **Fast** is better than *slow*:
   - **Vector** is better than *loop*
   - **Matrix** is better than *vector*
   - **Strided** is better than *scattered*
   - **Contiguous** is better than *strided*
   - **Broadcasting** is a great idea, use where possible
4. **Slow** is better than *unmaintainable*
5. Make it look like the **math**
6. **Special cases** aren't special enough to break rules...
7. Although **practicality** beats *purity*
8. **Pure** procedure is better than *impure*...
9. Although **practicality** beats *purity* again
10. **Private** is better than *public*
11. **Errors** should never pass *silently*...
12. Unless **errors** are explicitly *silenced*

## Standard Compliance

Standard-compliant code is **portable**:
- Free from particular OS/hardware architecture (ensures long-time code life, easy enabling on new hardware)
- Free from particular compiler vendor (can use different compilers for cross-checking/debugging)
- Use compiler flags to check standard adherence: `-std=f2018` (GNU), `-std18` (Intel)

## Readability Guidelines

**Use free source form**:
- 132 characters are better than 72 (use modern wide screens)
- For mathematical formulas, avoid splitting on multiple lines when possible
- Balance: don't fill each line to 132 chars, but don't limit yourself to 72

**Indentation**:
- Use consistent number of spaces (project standard: 3 spaces)
- Increase indentation when data scope changes
- Indent blocks within all control constructs
- Indent code after named constructs so names stand out
- Break lines in logical places; indent continued lines double the block indent
- Use blank lines to separate related parts of a program

**White spaces**:
- Avoid hard tabs
- Use spaces around operators: `if (foo > bar) then`
- Align similar code statements (declarations, comments)
- Place a space after all commas: `x(i, j, k) = foo(i, j, k:k+1)`

**Example of good style**:
```fortran
module well_formatted_module
implicit none
private
public :: foo_type

type :: foo_type
   private
   logical                       :: is_good = .false.
   character(len=:), allocatable :: name
   contains
      private
      procedure, public, pass(self) :: init
endtype foo_type

contains
   pure subroutine init(self, name)
   class(foo_type),        intent(inout) :: self
   character(*), optional, intent(in)    :: name

   if (allocated(self%name)) deallocate(self%name)
   check_name: if (present(name)) then
                  if (len_trim(adjustl(name)) == 0) exit check_name
                  self%name = trim(adjustl(name))
               else
                  self%name = 'J. Doe'
   endif check_name
   if (allocated(self%name)) self%is_good = .true.
   endsubroutine init
endmodule well_formatted_module
```

## Language Features and Standards

- Use modern Fortran (2003/2008/2018): `submodule`, `associate`, `block`, `error stop`
- Explicit interfaces via `module` and `submodule` (mandatory for all public procedures)
- `implicit none` everywhere (no exceptions)
- Use `iso_fortran_env` and `iso_c_binding` for portability
- All module entities should be `private` by default; explicitly declare `public` only what is part of the API
- Type components should be `private` unless external access is necessary

## Error Handling

- Never ignore error codes from allocate, MPI, I/O operations
- Use `stat=` and `errmsg=` with allocate/deallocate
- Check MPI return codes explicitly
- Provide meaningful error messages with context
- Use `error stop` with informative messages for unrecoverable errors

## Array Programming

- Prefer whole-array syntax when vectorizable: `a(:,:,:) = b(:,:,:) + c(:,:,:)`
- Use assumed-shape arrays for flexibility: `real(rp), intent(in) :: a(:,:,:)`
- Add `contiguous` attribute for performance-critical paths: `real(rp), intent(in), contiguous :: a(:,:,:)`
- Avoid assumed-size arrays: `real(rp) :: a(*)` (use assumed-shape or explicit-shape)
- Inner loops on first index for cache efficiency (column-major order)

## Type-Bound Procedures

- Use type-bound procedures for polymorphism and encapsulation
- Mark `procedure :: method => implementation` for clarity
- Use `nopass` only when justified (static-like methods)

## Preprocessor Usage

- Minimize preprocessor macros - prefer Fortran abstractions when possible
- When unavoidable, use for: compiler detection, backend selection, MPI/GPU conditionals
- Never use macros for code that can be written in standard Fortran

---

## Modern Fortran Patterns and Examples

*Practical examples extracted from the "ZEN of Fortran" presentation - use these patterns as reference when writing or reviewing code.*

### Variables and Constants

**Use named constants instead of magic numbers**:
```fortran
integer, parameter                           :: rows_number=5
integer, parameter                           :: cols_number=3
integer, dimension(rows_number, cols_number) :: array
```

**Always use `parameter` for constants**:
```fortran
integer,       parameter :: cars_number=100
real,          parameter :: pi=3.14159265
character(12), parameter :: greetings='Hello World!'
```

**Initialize in executable section to avoid implicit SAVE**:
```fortran
! WARNING: Initialization at declaration has implicit SAVE - value persists between calls!
real function kinetic_energy(v)
   real, dimension(:), intent(in) :: v
   integer                        :: i
   real                           :: ke

   ke = 0.0  ! Initialize HERE, in executable section
   do i = 1, size(v)
      ke = ke + v(i)**2
   enddo
   kinetic_energy = 0.5*ke
end function kinetic_energy
```

### Kind Specifications (Portability)

**Using `iso_fortran_env` (preferred)**:
```fortran
use, intrinsic :: iso_fortran_env, only : int32, real64

integer(kind=int32) :: a_good_parametric_integer
real(kind=real64)   :: a_good_parametric_real
```

**Using `selected_*_kind` (alternative)**:
```fortran
integer, parameter :: I4 = selected_int_kind(9)        ! 4-bytes
integer, parameter :: R8 = selected_real_kind(15,307)  ! 8-bytes

integer(kind=I4) :: a_good_parametric_integer
real(kind=R8)    :: a_good_parametric_real
```

**Always specify kind suffix for literal constants**:
```fortran
use, intrinsic :: iso_fortran_env, only : R8 => real64

! CORRECT: Kind-specified literals
real(kind=R8), parameter :: pi_greek = atan(1._R8)*4._R8

! WRONG: Default precision literals lose precision
real(kind=R8), parameter :: pi_bad = atan(1.0)*4.0  ! Computed in single precision!
```

### Array Best Practices

**Column-major access pattern (critical for performance)**:
```fortran
integer, parameter                           :: rows_number=5
integer, parameter                           :: cols_number=3
integer, dimension(rows_number, cols_number) :: array
integer                                      :: c, r

! CORRECT: Inner loop on first index (contiguous memory access)
do c=1, cols_number
   do r=1, rows_number
      print*, array(r, c)
   enddo
enddo
```

**Array constructors and elemental operations**:
```fortran
integer, parameter              :: elems_number=10
real, parameter                 :: delta=1.0/elems_number
real, dimension(0:elems_number) :: abscissa
real, dimension(0:elems_number) :: temperature
integer                         :: e

! Array constructor with implied do
abscissa = [(delta*e, e=0, elems_number)]

! Elemental function operates on entire array - vectorizable
temperature = sin(abscissa)/2.0
```

**Dynamic allocation**:
```fortran
real, allocatable, dimension(:) :: abscissa
real, allocatable, dimension(:) :: temperature
integer                         :: elems_number, e, alloc_stat
character(256)                  :: alloc_msg
real                            :: delta

elems_number = 100
delta = 1.0/elems_number

allocate(abscissa(0:elems_number), stat=alloc_stat, errmsg=alloc_msg)
if (alloc_stat /= 0) error stop 'Allocation failed: '//trim(alloc_msg)

allocate(temperature(0:elems_number), stat=alloc_stat, errmsg=alloc_msg)
if (alloc_stat /= 0) error stop 'Allocation failed: '//trim(alloc_msg)

abscissa = [(delta*e, e=0, elems_number)]
temperature = sin(abscissa)/2.0
```

### Subprogram Best Practices

**Keep It Simple (KISS) - max 50-100 lines per procedure**:
```fortran
subroutine do_something(a, b, c)
   ! ... implementation
end subroutine do_something

! Use keyword arguments for clarity
call do_something(a=a, b=b, c=c)
call do_something(a=d, b=e, c=f)
```

**Always specify intent**:
```fortran
subroutine compute(input, output, both)
   real, intent(in)    :: input   ! Read-only
   real, intent(out)   :: output  ! Write-only (undefined on entry)
   real, intent(inout) :: both    ! Read-write
end subroutine compute
```

**With `intent(out)` and derived types, assign ALL components**:
```fortran
subroutine assign(this)
   type(mytype), intent(out) :: this
   this%x = 2
   this%y = 2.0  ! Assign ALL components - intent(out) makes them undefined!
end subroutine assign
```

**Handle optional arguments safely**:
```fortran
subroutine print_char(this, header)
   character(len=*),  intent(in) :: this
   logical, optional, intent(in) :: header

   if (present(header)) then
      if (header) print *, 'This is the header'
   endif
   print *, this
end subroutine print_char
```

### Module Best Practices

**Template for well-structured modules**:
```fortran
module library
   use another_library, only : external_concern => what_really_concern  ! Rename to avoid collisions
   implicit none
   private              ! Hide everything by default
   public :: foo        ! Expose only what's needed
   protected :: read_only_data  ! Read-only from outside module
contains
   pure subroutine foo(arg1, arg2, arg3)
      real, intent(in)  :: arg1
      real, intent(in)  :: arg2(:)
      real, intent(out) :: arg3
      ! Implementation
   end subroutine foo
end module library
```

### I/O Best Practices

**Use standard units from `iso_fortran_env`**:
```fortran
use, intrinsic :: iso_fortran_env, only : error_unit, input_unit, output_unit

read(input_unit, *) user_input
if (error_occur) then
   write(error_unit, '(A)') error_message
else
   write(output_unit, '(A)') output_message
endif
```

### Generic Interfaces

**Apply same procedure to different types**:
```fortran
module generic
   implicit none
   private
   public :: string

   interface string
      module procedure string_integer, string_real
   end interface
contains
   elemental function string_integer(n) result(str)
      integer, intent(in) :: n
      character(11)       :: str
      write(str, '(I11)') n
   end function string_integer

   elemental function string_real(n) result(str)
      real, intent(in) :: n
      character(13)    :: str
      write(str, '(E13.6E2)') n
   end function string_real
end module generic

! Usage:
! print*, string(42)      ! Calls string_integer
! print*, string(3.14)    ! Calls string_real
```

### Operator Overloading

**Create custom algebra for derived types**:
```fortran
module vector_t
   implicit none
   private
   public :: vector, operator(+), operator(.cross.)

   type :: vector
      real :: x = 0.0
      real :: y = 0.0
      real :: z = 0.0
   end type vector

   interface operator(+)
      module procedure add
   end interface

   interface operator(.cross.)
      module procedure crossproduct
   end interface
contains
   elemental function add(left, right) result(summ)
      type(vector), intent(in) :: left, right
      type(vector)             :: summ
      summ%x = left%x + right%x
      summ%y = left%y + right%y
      summ%z = left%z + right%z
   end function add

   elemental function crossproduct(vec1, vec2) result(cross)
      type(vector), intent(in) :: vec1, vec2
      type(vector)             :: cross
      cross%x = (vec1%y * vec2%z) - (vec1%z * vec2%y)
      cross%y = (vec1%z * vec2%x) - (vec1%x * vec2%z)
      cross%z = (vec1%x * vec2%y) - (vec1%y * vec2%x)
   end function crossproduct
end module vector_t

! Usage:
! type(vector) :: v1, v2, v3
! v3 = v1 + v2
! v3 = v1 .cross. v2
```

### Object-Oriented Programming

**Encapsulation - hide implementation, expose interface**:
```fortran
module gas_t
   use vecfor, only : vector
   implicit none
   private
   public :: gas

   type :: gas
      private  ! Hide all members
      real         :: Cp_ = 0.0
      real         :: Cv_ = 0.0
      real         :: density_ = 0.0
      type(vector) :: momentum_
   contains
      private
      procedure, public :: initialize
      procedure, public :: density       ! Getter
      procedure, public :: temperature   ! Computed property
      procedure, public :: speed_of_sound
   end type gas
contains
   ! Implementation...
end module gas_t

! Usage:
! use gas_t, only : gas
! type(gas) :: air
! call air%initialize(Cp=1004.8, Cv=716.0, density=0.125)
! if (air%speed() > air%speed_of_sound()) then
!    ! Supersonic flow
! endif
```

**Inheritance - extend existing classes**:
```fortran
module gas_multifluid_t
   use gas_t, only : gas
   implicit none
   private
   public :: gas_multifluid

   type, extends(gas) :: gas_multifluid
      private
      integer           :: species_number = 0
      real, allocatable :: concentrations(:)
   contains
      procedure, public :: initialize  ! Override parent
   end type gas_multifluid
end module gas_multifluid_t
```

**Polymorphism with abstract types**:
```fortran
type, abstract :: integrand
contains
   procedure(time_derivative), deferred, public :: t
   generic, public :: operator(+) => add
   generic, public :: operator(*) => multiply
   generic, public :: assignment(=) => assign_integrand
end type integrand

! Solvers work with ANY extension of integrand:
subroutine integrate(U, Dt, t)
   class(integrand), intent(inout) :: U   ! Polymorphic!
   real,             intent(in)    :: Dt
   real, optional,   intent(in)    :: t

   U = U + U%t(t=t) * Dt  ! Works for any integrand extension
end subroutine integrate
```

### OpenMP Parallelization Patterns

**Basic parallel loop**:
```fortran
use omp_lib

!$omp parallel shared(a, b, c) private(i)
!$omp do
do i = 1, n
   c(i) = a(i) + b(i)
end do
!$omp end do
!$omp end parallel
```

**Use `default(none)` for safety - forces explicit scoping**:
```fortran
!$omp parallel default(none) shared(a, b, c, n) private(i)
! All variables must be explicitly scoped - catches bugs
!$omp end parallel
```

**Reduction for accumulation**:
```fortran
sum = 0
!$omp parallel do reduction(+:sum)
do i = 1, n
   sum = sum + x(i)
end do
!$omp end parallel do
```

**Use `workshare` for array syntax**:
```fortran
!$omp parallel shared(aa, bb, cc, dd)
!$omp workshare
cc = aa * bb
dd = aa + bb
!$omp end workshare
!$omp end parallel
```

## Quick Reference Table

| Topic | Best Practice |
|-------|---------------|
| Typing | Always use `implicit none` |
| Constants | Use `parameter` attribute, avoid magic numbers |
| Array access | Iterate inner index fastest (column-major) |
| Kinds | Use `iso_fortran_env` or `selected_*_kind` |
| Literal constants | Add kind suffix: `3.14_R8` |
| Initialization | Initialize in executable section (avoid implicit SAVE) |
| Optional arguments | Check `present()` before accessing |
| `intent(out)` | Assign all derived type components |
| Modules | Use `private` by default, expose with `public` |
| OpenMP | Use `default(none)`, `reduction` for accumulation |

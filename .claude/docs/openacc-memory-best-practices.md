# OpenACC Memory Management Best Practices

This guide covers memory management patterns for OpenACC GPU programming in Fortran, with specific focus on avoiding common pitfalls.

## Variable Types and Device Handling

| Type | Example | Needs `!$acc declare`? | Notes |
|------|---------|------------------------|-------|
| Scalar `parameter` | `integer, parameter :: N = 5` | **No** | Value inlined at compile time |
| Array `parameter` | `real, parameter :: C(5,5) = ...` | **Yes** | `!$acc declare copyin(C)` |
| Module `allocatable` | `real, allocatable :: arr(:)` | **Yes** | `!$acc declare create(arr)` |
| Module fixed-size | `real :: arr(100)` | **Yes** | `!$acc declare create(arr)` |
| Local variable | `real :: arr(100)` in subroutine | **No** | Use `private` clause |
| Subroutine argument | `real, intent(in) :: arr(:)` | **No** | Use `deviceptr` or data clauses |

## Module-Level Variables

### Scalar Parameters (No action needed)
```fortran
module my_module
   integer, parameter :: S_MAX = 5  ! Inlined - no device handling needed
end module
```

### Array Parameters (Need copyin declaration)
```fortran
module coefficients
   use penf
   implicit none

   real(R8P), parameter :: FD1_CC(5,5) = reshape([...], [5,5])
   !$acc declare copyin(FD1_CC)  ! Make available on device
end module
```

### Allocatable Module Arrays
```fortran
module data_module
   use penf
   implicit none

   real(R8P), allocatable :: global_array(:,:,:)
   !$acc declare create(global_array)  ! Declare device copy exists

contains
   subroutine initialize(n)
      integer, intent(in) :: n
      allocate(global_array(n,n,n))
      global_array = 0.0_R8P
      !$acc update device(global_array)  ! Copy data to device
   end subroutine

   subroutine cleanup()
      !$acc exit data delete(global_array)
      deallocate(global_array)
   end subroutine
end module
```

## Local Variables and Private Clause

### Correct: Local variables can be private
```fortran
subroutine compute(q_gpu, n, s)
   real(R8P), intent(inout) :: q_gpu(:,:,:)
   integer, intent(in) :: n, s
   real(R8P) :: stencil(1-s:1+s)  ! LOCAL variable
   real(R8P) :: result
   integer :: i, j, k

   !$acc parallel loop collapse(3) private(stencil, result) deviceptr(q_gpu)
   do k = 1, n
   do j = 1, n
   do i = 1, n
      stencil(:) = q_gpu(i-s:i+s, j, k)
      ! ... compute with stencil
   end do
   end do
   end do
end subroutine
```

### Wrong: Module-level variables cannot be private
```fortran
module bad_example
   real(R8P) :: stencil(10)  ! MODULE-level - ERROR!
contains
   subroutine compute()
      !$acc parallel loop private(stencil)  ! FAILS: No device symbol
      ! ...
   end subroutine
end module
```

**Error message:** `No device symbol for address reference`

## Array Sections in Device Code

### Problem: Non-contiguous array sections as arguments

Passing non-contiguous array sections to device subroutines can cause illegal memory access:

```fortran
! DANGEROUS: Non-contiguous array section passed to device routine
!$acc parallel loop deviceptr(q_gpu)
do i = 1, n
   call device_routine(q_gpu(i, 1:m, k, :))  ! May fail!
end do
```

### Solution 1: Copy to local contiguous array
```fortran
!$acc parallel loop private(local_buffer) deviceptr(q_gpu)
do i = 1, n
   local_buffer(:) = q_gpu(i, 1:m, k, var)  ! Copy to contiguous local
   call device_routine(local_buffer)
end do
```

### Solution 2: Pass full array with scalar indices
```fortran
pure subroutine device_routine(arr, i, j, k, n, ngc)
   real(R8P), intent(in) :: arr(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer, intent(in) :: i, j, k, n, ngc
   integer :: m
   !$acc routine seq

   ! Access directly with scalar indices
   do m = 1, n
      ! Use arr(b, i+m, j, k, var) directly
   end do
end subroutine
```

## Memory Layout Best Practices

### Fortran Column-Major Order

Fortran stores arrays in column-major order. The leftmost index varies fastest in memory:

```fortran
! Memory layout: arr(1,1), arr(2,1), arr(3,1), arr(1,2), arr(2,2), ...
real :: arr(ni, nj, nk)
```

### Optimal Loop Ordering for CPU
```fortran
! GOOD: Innermost loop on leftmost index (stride-1 access)
do k = 1, nk
   do j = 1, nj
      do i = 1, ni  ! Fastest varying - contiguous in memory
         arr(i, j, k) = ...
      end do
   end do
end do
```

### GPU Coalesced Access

For GPU, adjacent threads should access adjacent memory locations:

```fortran
! GOOD for GPU: Collapse outer loops, let threads handle inner indices
!$acc parallel loop collapse(3) gang vector
do k = 1, nk
do j = 1, nj
do i = 1, ni
   ! Thread (i,j,k) accesses arr(i,j,k)
   ! Adjacent threads access adjacent i values = coalesced
   arr(i, j, k) = ...
end do
end do
end do
```

### ADAM Field Array Convention

ADAM uses 5D arrays: `(nb, ni, nj, nk, nv)` where:
- `nb` = block index (first dimension for efficient block-parallel access)
- `ni, nj, nk` = spatial indices (with ghost cells: `1-ngc:n+ngc`)
- `nv` = variable index

```fortran
real(R8P) :: q_gpu(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv)
```

## Data Movement

### Explicit Data Regions
```fortran
!$acc data copyin(input) copyout(output) create(temp)
   !$acc parallel loop
   do i = 1, n
      temp(i) = input(i) * 2.0
   end do

   !$acc parallel loop
   do i = 1, n
      output(i) = temp(i) + 1.0
   end do
!$acc end data
```

### Unstructured Data Lifetime
```fortran
subroutine init()
   allocate(arr(n))
   !$acc enter data create(arr)
end subroutine

subroutine compute()
   !$acc parallel loop present(arr)
   do i = 1, n
      arr(i) = ...
   end do
end subroutine

subroutine finalize()
   !$acc exit data delete(arr)
   deallocate(arr)
end subroutine
```

### Device Pointers (FUNDAL pattern)
```fortran
! Array allocated directly on device via acc_malloc / dev_alloc
real(R8P), pointer :: q_gpu(:,:,:,:,:)
call dev_alloc(fptr_dev=q_gpu, ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1])

! Use deviceptr to tell compiler it's already on device
!$acc parallel loop deviceptr(q_gpu)
do i = 1, n
   q_gpu(b,i,j,k,v) = ...
end do
```

## Common Pitfalls

### 1. Assumed-shape arrays with non-default lower bounds
```fortran
! May have issues in device code
subroutine problematic(arr)
   real(R8P), intent(in) :: arr(1-s:)  ! Lower bound from dummy argument
   !$acc routine seq
```

**Solution:** Pass bounds explicitly and use them for indexing.

### 2. Associate variables in OpenACC regions
```fortran
! May not work with firstprivate
associate(x => self%data%x)
   !$acc parallel loop firstprivate(x)  ! Can fail with some compilers
```

**Solution:** Use `copyin` instead of `firstprivate` for associate variables.

### 3. Forgetting to update device after host modification
```fortran
arr(:) = new_values  ! Modified on host
! WRONG: Device copy is stale
!$acc parallel loop present(arr)

! CORRECT: Update device first
!$acc update device(arr)
!$acc parallel loop present(arr)
```

### 4. Pointer vs allocatable arrays
```fortran
! Allocatable - compiler manages device allocation
real(R8P), allocatable :: arr(:)
!$acc declare create(arr)

! Pointer - YOU manage device memory (use FUNDAL dev_alloc)
real(R8P), pointer :: arr_gpu(:)
! Must use deviceptr, not present/create
```

## Debugging Tips

1. **Enable runtime checks:** `NV_ACC_DEBUG=1` environment variable
2. **Check data presence:** `acc_is_present(arr, size)`
3. **Synchronize for error detection:** `!$acc wait` forces synchronization and surfaces async errors
4. **Use compute-sanitizer:** `compute-sanitizer ./executable` catches memory errors

## References

- [OpenACC Specification](https://www.openacc.org/specification)
- [NVIDIA HPC SDK Documentation](https://docs.nvidia.com/hpc-sdk/)
- [OpenACC Best Practices Guide](https://www.openacc.org/sites/default/files/inline-files/OpenACC_Programming_Guide_0.pdf)

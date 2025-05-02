module adam_prism_physics_object

    use penf, only : I4P, R8P

    implicit none

    integer(I4P), parameter :: VAR_DX = 1_I4P
    integer(I4P), parameter :: VAR_DY = 2_I4P
    integer(I4P), parameter :: VAR_DZ = 3_I4P
    integer(I4P), parameter :: VAR_BX = 4_I4P
    integer(I4P), parameter :: VAR_BY = 5_I4P
    integer(I4P), parameter :: VAR_BZ = 6_I4P
    integer(I4P), parameter :: VAR_JX = 7_I4P
    integer(I4P), parameter :: VAR_JY = 8_I4P
    integer(I4P), parameter :: VAR_JZ = 9_I4P
    type :: prism_physics_object
    !< PRISM physics class definition.
    type(mpih_object)                   :: mpih         !< MPI handler.
    integer(I4P)                        :: nv=9_I4P     !< Number of variables (see below).
    !integer(I4P)                        :: ns=1_I4P     !< Number of species.
    !integer(I4P)                        :: nv_aux=9_I4P !< Number of auxiliary variables (rns+r+u+v+w+p+g=ns+6).
    !integer(I4P)                        :: np=5_I4P     !< Number of 1D primitive variables (rns+r+un+p+g=ns+4).
    !type(nasto_eos_object), allocatable :: eos(:)       !< Equations of state of each specie [1:ns].
    !real(R8P), allocatable :: q_field(:,:,:,:,:) 

    ! conservative variables are arranged as follows (and they coincide with auxiliary ones):
    ! q_field(1,i,j,k,b): Dx
    ! q_field(2,i,j,k,b): Dy
    ! q_field(3,i,j,k,b): Dz
    ! q_field(4,i,j,k,b): Bx
    ! q_field(5,i,j,k,b): By
    ! q_field(6,i,j,k,b): Bz
!    contains
!       ! public methods
!       !procedure, pass(self) :: description            !< Return pretty-printed object description.
!       procedure, pass(self) :: initialize             !< Initialize physics.
!       !procedure, pass(self) :: load_from_file         !< Load config from file.
!       !procedure, pass(self) :: conservative2primitive !< Return primitive variables from conservative ones.
!       !procedure, pass(self) :: primitive2conservative !< Return conservative variables from primitive ones.
    endtype prism_physics_object
! contains
! subroutine initialize(self, ni, nj, nk, nb)
!    class(prism_physics_object) :: self
!    integer(I4P), intent(in) :: ni, nj, nk, nb
!    allocate(self%q_field(nv,ni,nj,nk,nb))
!    self%q_field = 0._R8P
! endsubroutine

endmodule adam_prism_phisics_object

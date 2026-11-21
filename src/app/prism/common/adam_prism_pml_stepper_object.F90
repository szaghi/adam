!< PRISM, the PML auxiliary fields as SELAH companions of the field stepper (#134 D5).

#include "fundal.H"

module adam_prism_pml_stepper_object
!< PRISM, the PML auxiliary fields as SELAH companions of the field stepper (#134 D5).
!<
!< The auxiliary fields psi of a perfectly matched layer live in face-local arrays, one per active face of the domain
!< (`q_pml_x_m`, ..., `q_pml_z_p`: `[4, cells, nj, nk, blocks]` on the CPU, `[blocks, cells, nj, nk, 4]` on the FNL
!< backend). They advance in the same step as the fields, their residual built inside the fields' (psi' = g - lambda
!< psi, g from the field derivatives): each face is a SELAH companion of the field stepper (`initialize_companion`:
!< the same scheme, form and stages on q), driven through the same stage window as the field (`begin_step`,
!< `begin_stage(k)`, the residual into `rhs_register(f)`, `end_stage(k)`, `finish_step`). SELAH owns the stage states,
!< the registers and the update; PRISM writes the residual. Stepped together, field and companion are the two stepped
!< as one state (szaghi/SELAH#4).
!<
!< **Exact decay** (`[PML] exact_decay`, default false; szaghi/SELAH#4 variant F): with it the residual is the source g
!< alone and lambda = alpha + gamma / kappa (written by the residual kernels into `decay(f)`, static) is integrated
!< exactly by SELAH (`decay=` at every begin_stage and at finish_step: exponential Euler on the Runge-Kutta-averaged
!< source), stable for any lambda dt; without it the explicit step of psi' = g - lambda psi, bitwise the former
!< integrator's arithmetic but for the order of the stage sums. The first stage of an explicit scheme has c = 0, so
!< the decay of stage 1 is never read before the first residual writes it.
!<
!< The same source serves the CPU (host arrays) and the FNL backend (device arrays, SELAH's device kernels): the face
!< arrays are passed by the caller, the registers and `decay(f)` are allocated through FUNDAL.
use :: adam_mpih_global,    only : mpih
use :: adam_stepper_object, only : stepper_object
use :: fundal,              only : dev_alloc, dev_free
use :: penf,                only : I4P, R8P
use :: selah,               only : stepper, box, box_interior

implicit none
private
public :: prism_pml_stepper_object

integer(I4P), parameter :: FACES = 6_I4P !< Faces of the domain.

type :: face_decay
   !< Decay rate of one face (exact decay only).
   real(R8P), pointer, contiguous :: lam(:,:,:,:,:) => null() !< lambda = alpha + gamma / kappa, the face's shape.
endtype face_decay

type :: prism_pml_stepper_object
   !< The six faces' companions.
   logical           :: enabled = .false.     !< A PML integrated with the field stages.
   logical           :: exact_decay = .false. !< The decay integrated exactly (SELAH `decay=`).
   logical           :: active(FACES) = .false. !< The face has a layer.
   type(stepper)     :: st(FACES)             !< Companions.
   type(box)         :: bx(FACES)             !< Boxes: every element of the face array.
   type(face_decay)  :: decay(FACES)          !< Decay rates (exact decay).
   contains
      procedure, pass(self) :: begin_stage   !< Open stage k of face f.
      procedure, pass(self) :: begin_step    !< Open the step of face f.
      procedure, pass(self) :: destroy       !< Free the companions.
      procedure, pass(self) :: end_stage     !< Close stage k of face f.
      procedure, pass(self) :: finish_step   !< Close the step of face f.
      procedure, pass(self) :: initialize    !< Initialize the companions of the active faces.
      procedure, pass(self) :: rhs_register  !< Residual register of face f (stage window).
      procedure, pass(self) :: stage_state   !< State of the active stage of face f (stage window).
      procedure, pass(self) :: zero_rhs      !< Zero the residual register of face f (stage window).
endtype prism_pml_stepper_object

contains
   subroutine initialize(self, ts, shapes, active, exact_decay)
   !< Initialize a companion of the field stepper per active face, of the face array's shape.
   class(prism_pml_stepper_object), intent(inout) :: self              !< PML companions.
   type(stepper_object),            intent(in)    :: ts                !< The field's stepper (initialized).
   integer(I4P),                    intent(in)    :: shapes(5,FACES)   !< Shape of each face array.
   logical,                         intent(in)    :: active(FACES)     !< The face has a layer.
   logical,                         intent(in)    :: exact_decay       !< Integrate the decay exactly.
   integer(I4P)                                   :: f, ierr           !< Counter, error status.

   call self%destroy
   if (.not.ts%is_butcher()) call mpih%error_stop(msg=': the PML auxiliary fields advance with the field stages of '// &
                                                      'a Runge-Kutta scheme of the Butcher form (runge-kutta-ssp-*)')
   self%enabled = .true.
   self%exact_decay = exact_decay
   self%active = active
   do f=1, FACES
      if (.not.self%active(f)) cycle
      call self%st(f)%initialize_companion(ts%st, lb=[1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P], ub=shapes(:,f))
      ! every element evolves (the elements beyond a block's cells have a zero residual and stay zero)
      self%bx(f) = box_interior(lb=[1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P], ub=shapes(:,f), ngc=0_I4P, components=.false.)
      if (self%exact_decay) then
         call dev_alloc(fptr_dev=self%decay(f)%lam, lbounds=[1_I4P, 1_I4P, 1_I4P, 1_I4P, 1_I4P], ubounds=shapes(:,f), &
                        init_value=0._R8P, ierr=ierr)
         if (ierr /= 0) call mpih%error_stop(msg=': prism_pml_stepper_object: the decay rates could not be allocated')
      endif
   enddo
   endsubroutine initialize

   subroutine destroy(self)
   !< Free the companions and the decay rates.
   class(prism_pml_stepper_object), intent(inout) :: self !< PML companions.
   integer(I4P)                                   :: f    !< Counter.

   do f=1, FACES
      call self%st(f)%destroy
      if (associated(self%decay(f)%lam)) call dev_free(self%decay(f)%lam)
      nullify(self%decay(f)%lam)
   enddo
   self%enabled = .false. ; self%exact_decay = .false. ; self%active = .false.
   endsubroutine destroy

   subroutine begin_step(self, f, q)
   !< Open the step of face f.
   class(prism_pml_stepper_object), intent(inout) :: self          !< PML companions.
   integer(I4P),                    intent(in)    :: f             !< Face.
   real(R8P),                       intent(in)    :: q(:,:,:,:,:)  !< Face array.

   call self%st(f)%begin_step(self%bx(f), q)
   endsubroutine begin_step

   subroutine begin_stage(self, f, k, dt, q)
   !< Open stage k of face f: SELAH builds its state (with the exact decay, from the decay rates).
   class(prism_pml_stepper_object), intent(inout) :: self          !< PML companions.
   integer(I4P),                    intent(in)    :: f             !< Face.
   integer(I4P),                    intent(in)    :: k             !< Stage.
   real(R8P),                       intent(in)    :: dt            !< Time step.
   real(R8P),                       intent(in)    :: q(:,:,:,:,:)  !< Face array.

   if (self%exact_decay) then
      call self%st(f)%begin_stage(k, self%bx(f), dt, q, decay=self%decay(f)%lam)
   else
      call self%st(f)%begin_stage(k, self%bx(f), dt, q)
   endif
   endsubroutine begin_stage

   subroutine end_stage(self, f, k, dt, q)
   !< Close stage k of face f, its residual written into `rhs_register(f)`.
   class(prism_pml_stepper_object), intent(inout) :: self          !< PML companions.
   integer(I4P),                    intent(in)    :: f             !< Face.
   integer(I4P),                    intent(in)    :: k             !< Stage.
   real(R8P),                       intent(in)    :: dt            !< Time step.
   real(R8P),                       intent(inout) :: q(:,:,:,:,:)  !< Face array.

   call self%st(f)%end_stage(k, self%bx(f), dt, q)
   endsubroutine end_stage

   subroutine finish_step(self, f, dt, q)
   !< Close the step of face f: SELAH advances the face array (with the exact decay, from the decay rates).
   class(prism_pml_stepper_object), intent(inout) :: self          !< PML companions.
   integer(I4P),                    intent(in)    :: f             !< Face.
   real(R8P),                       intent(in)    :: dt            !< Time step.
   real(R8P),                       intent(inout) :: q(:,:,:,:,:)  !< Face array.

   if (self%exact_decay) then
      call self%st(f)%finish_step(self%bx(f), dt, q, decay=self%decay(f)%lam)
   else
      call self%st(f)%finish_step(self%bx(f), dt, q)
   endif
   endsubroutine finish_step

   function rhs_register(self, f) result(r)
   !< Return the residual register of face f (null outside its stage window).
   class(prism_pml_stepper_object), intent(in) :: self          !< PML companions.
   integer(I4P),                    intent(in) :: f             !< Face.
   real(R8P), pointer                          :: r(:,:,:,:,:)  !< Residual register.

   r => self%st(f)%rhs_register()
   endfunction rhs_register

   function stage_state(self, f) result(u)
   !< Return the state of the active stage of face f (null outside its stage window).
   class(prism_pml_stepper_object), intent(in) :: self          !< PML companions.
   integer(I4P),                    intent(in) :: f             !< Face.
   real(R8P), pointer                          :: u(:,:,:,:,:)  !< Stage state.

   u => self%st(f)%stage_state()
   endfunction stage_state

   subroutine zero_rhs(self, f)
   !< Zero the residual register of face f: the residual kernels write the elements of each block's cells only, the
   !< others keep a zero residual (and a zero state).
   class(prism_pml_stepper_object), intent(in) :: self          !< PML companions.
   integer(I4P),                    intent(in) :: f             !< Face.
   real(R8P), pointer                          :: r(:,:,:,:,:)  !< Residual register.
   integer(I4P)                                :: i1, i2, i3, i4, i5 !< Counters.

   r => self%st(f)%rhs_register()
   !$acc parallel loop independent collapse(5) DEVICEVAR(r)
   !$omp OMPLOOP collapse(5) DEVICEPTR(r)
   do i5=1, ubound(r, 5)
      do i4=1, ubound(r, 4)
         do i3=1, ubound(r, 3)
            do i2=1, ubound(r, 2)
               do i1=1, ubound(r, 1)
                  r(i1,i2,i3,i4,i5) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine zero_rhs
endmodule adam_prism_pml_stepper_object

!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, CPU backend.
module adam_prism_cpu_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, CPU backend.

! ADAM modules
use adam_common_library
! PRISM modules
use adam_prism_common_library
! third party modules
use penf
use mpi

implicit none
private
public :: prism_cpu_object

! pointer (abstract) procedures
procedure(compute_convective_fluxes_interface), pointer :: compute_fluxes_Maxwell=>null() !< Compute convective fluxes.

type, extends(prism_common_object) :: prism_cpu_object !commentate procedure AMR e IB
   !< Maxwell equations system class definition, CPU backend.
   real(R8P), allocatable :: flxyz_c(:,:,:,:,:,:,:) !< Fluxes at cell center with +/- decomposition for all directions.
   real(R8P), allocatable ::   flx_f(:,:,:,:,:    ) !< Fluxes along x at cell face.
   real(R8P), allocatable ::   fly_f(:,:,:,:,:    ) !< Fluxes along y at cell face.
   real(R8P), allocatable ::   flz_f(:,:,:,:,:    ) !< Fluxes along z at cell face.
   !< Pointer (abstract) TBP.
   procedure(compute_residuals_interface), pass(self),pointer :: compute_residuals=>null()!< Compute residuals, space operator.
   procedure(integrate_interface),         pass(self),pointer :: integrate        =>null()!< Integrate, time operator.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu !< Allocate CPU data.
      procedure, pass(self) :: initialize   !< Initialize the equation.
      ! IO methods
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC/sources
      procedure, pass(self) :: apply_fWL_correction              !< Apply fWLayer correction (if present)
      procedure, pass(self) :: compute_coils_current             !< Compute current coils sources.
      procedure, pass(self) :: impose_div_coil_correction        !< Impose coil divergence correction.
      procedure, pass(self) :: set_boundary_conditions           !< Set boundary conditions of equation.
      procedure, pass(self) :: compute_residuals_BC
      procedure, pass(self) :: update_q_BC
      procedure, pass(self) :: set_initial_conditions            !< Set initial conditions (and coils) of equation.
      procedure, pass(self) :: update_ghost                      !< Update ghost cells and set boundary conditions.
      procedure, pass(self) :: compute_field_mean_value          !< Compute field mean value.
      ! numerical methods
      procedure, pass(self) :: compute_auxiliary_fields !< Compute auxiliary fields.
      procedure, pass(self) :: compute_dt               !< Compute time step.
      procedure, pass(self) :: compute_energy           !< Compute energy.
      procedure, pass(self) :: compute_energy_error     !< Compute energy error.
      procedure, pass(self) :: impose_div_free          !< Impose divergence-free property.
      procedure, pass(self) :: impose_ct_correction     !< Impose Constrained Transport correction on q(ivar:ivar+2).
      procedure, pass(self) :: simulate                 !< Perform the simulation.
endtype prism_cpu_object

interface
   subroutine compute_residuals_interface(self, q, dq, s)
   !< Compute residuals of equation, space operator.
   import :: prism_cpu_object, R8P, I4P
   class(prism_cpu_object), intent(inout) :: self   !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)  !< Conservative variables.
   real(R8P),               intent(inout) :: dq(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:) !< Residuals.
   integer(I4P),  optional, intent(in)    :: s !< Stage counter.
   endsubroutine compute_residuals_interface

   subroutine integrate_interface(self)
   !< Integrate equation, time operator.
   import :: prism_cpu_object, R8P
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_interface
endinterface

contains
   ! auxiliary methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.

   call self%mpih%print_message('prism_cpu_object%allocate_cpu start')
   msg_ = self%mpih%myrankstr//'prism_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   msg = msg_//' flxyz_c '
   call allocate_variable(var=self%flxyz_c,ulb=reshape([1,nv, &
                                                        1,3,  & ! f at center, f positive at center, f negative at center
                                                        1,3,  & ! direction x, y, z
                                                        1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,7]),msg=msg)
   self%flxyz_c = 0._R8P
   msg = msg_//' flx_f '
   call allocate_variable(var=self%flx_f,ulb=reshape([1,nv,0,ni,1,nj,1,nk,1,nb],[2,5]),msg=msg)
   self%flx_f = 0._R8P
   msg = msg_//' fly_f '
   call allocate_variable(var=self%fly_f,ulb=reshape([1,nv,1,ni,0,nj,1,nk,1,nb],[2,5]),msg=msg)
   self%fly_f = 0._R8P
   msg = msg_//' flz_f '
   call allocate_variable(var=self%flz_f,ulb=reshape([1,nv,1,ni,1,nj,0,nk,1,nb],[2,5]),msg=msg)
   self%flz_f = 0._R8P
   endassociate
   call self%mpih%print_message('prism_cpu_object%allocate_cpu finish')
   endsubroutine allocate_cpu

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(prism_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih%initialize(do_mpi_init=.true., verbose=.true.)
   call self%mpih%print_message('prism_cpu_object%initialize start')
   call self%initialize_common(field = self%adam%field, filename=filename, memory_avail=self%mpih%memory_avail)
   call self%allocate_cpu

   ! set pointer (abstract) TBP
   if (self%physics%physical_model == EM_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_BLANES_MOAN)
         self%integrate => integrate_blanesmoan
      case(NUM_SCHEME_TIME_CFM)
         self%integrate => integrate_cfm
      case(NUM_SCHEME_TIME_LEAPFROG)
         self%integrate => integrate_leapfrog
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%rk%scheme)
         case(RK_1, RK_2, RK_3)
            self%integrate => integrate_rk_ls
         case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
            self%integrate => integrate_rk_ssp
         case(RK_YOSHIDA)
            self%integrate => integrate_rk_yoshida
         endselect
      endselect
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then !Metterei qualche error stop sulle combinazioni non valide
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_LEAPFROG)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
            !self%integrate =>
         endselect
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
         select case(self%rk_pic%scheme)
         case(RK_1, RK_2, RK_3)
            !self%integrate => integrate_rk_ls_pic
         case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
            self%integrate => integrate_rk_ssp_pic
         case(RK_YOSHIDA)
            !self%integrate => integrate_rk_yoshida_pic
         endselect
         endselect
      endselect
   endif

   select case(self%numerics%scheme_space)
   case(NUM_SCHEME_SPACE_WENO)
      self%compute_residuals   => compute_residuals_weno
   case(NUM_SCHEME_SPACE_FD_CENTERED)
      self%compute_residuals   => compute_residuals_fd_centered
   case(NUM_SCHEME_SPACE_FV_CENTERED)
      self%compute_residuals   => compute_residuals_fv_centered
   endselect

   select case(self%numerics%div_corr_var)
   case(DIV_CORR_VAR_POISS)
      compute_fluxes_Maxwell => compute_convective_fluxes_maxwell
   case(DIV_CORR_VAR_HYPER)
      if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
         compute_fluxes_Maxwell => compute_convective_fluxes_maxwell_div_d
      elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
         compute_fluxes_Maxwell => compute_convective_fluxes_maxwell_div_b
      elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
         compute_fluxes_Maxwell => compute_convective_fluxes_maxwell_div_d_b
      endif
   case default
      compute_fluxes_Maxwell => compute_convective_fluxes_maxwell
   endselect

   print '(A)', self%mpih%description()
   call self%mpih%print_message('prism_cpu_object%initialize finish')
   endsubroutine initialize

   ! IO methods
   subroutine save_residuals(self)
   !< Save residuals history.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call self%field%compute_normL2_residuals(dq=self%dq, norm=self%field%residuals)
      do v=1, self%nv
         call MPI_ALLREDUCE(MPI_IN_PLACE, self%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
         self%field%residuals(v) = sqrt(self%field%residuals(v))/sqrt(real(self%ni*self%nj*self%nk, R8P))
      enddo
      if (self%mpih%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                           blocks_number=self%blocks_number, residuals=self%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q=self%q)
      call self%compute_auxiliary_fields

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max)) then
         call self%slices%save_mat(basename=self%io%output_basename, &
                                   it=self%time%it,                  &
                                   it_max=self%time%it_max,          &
                                   time=self%time%time,              &
                                   time_max=self%time%time_max,      &
                                   adam=self%adam,                   &
                                   q=self%q,                         &
                                   q_name=self%q_name)

      endif
   endif
   if (self%pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      call write_single_particle_output(filename='single_particle_output.dat', time=self%time%time, q_pic=self%q_pic)
   endif
   endsubroutine save_simulation_data

   ! IC/BC/sources
   !subroutine compute_coils_current(self, q, gamma)
   !!< Compute current coils sources.
   !class(prism_cpu_object), intent(inout)        :: self             !< The equation.
   !real(R8P),               intent(inout)        :: q(1:,          &
   !                                                   1-self%ngc:, &
   !                                                   1-self%ngc:, &
   !                                                   1-self%ngc:, &
   !                                                   1:)            !< Conservative variables.
   !real(R8P),               intent(in), optional :: gamma            !< RK coefficient.
   !real(R8P)                                     :: current_density  !< Current density.
   !real(R8P)                                     :: current_density_o!< Current density.
   !real(R8P)                                     :: g                !< Starting polynomial transitory of coils.
   !real(R8P)                                     :: time_s           !< Local time.
   !integer(I4P)                                  :: w_, w_c_         !< Step function coeff to avoid if in parallel regions.
   !real(R8P)                                     :: g_, f_           !< Current coefficients.
   !integer(I4P)                                  :: coil_id          !< Uniq coild ID.
   !integer(I4P)                                  :: i,j,k,b,n        !< Counter.
!
   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number,   &
   !          time=>self%time%time, A=>self%coil%A, f=>self%coil%f, phase=>self%coil%phase,              &
   !          coil_flag =>self%coil%coil_flag, td=>self%coil%td, J_vec=>self%coil%J_vec,                 &
   !          var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz,     &
   !          dx=>self%field%dxyz(1,1), dt=>self%time%dt)
!
   !if (present(gamma)) then
   !   time_s = time + dt*gamma
   !else
   !   time_s = time
   !end if
   !if (self%coil%total_coils_number >= 1_I4P) then
   !   ! Azzero termini sorgenti (con il PIC questa inizializzazione non va bene, bisogna ragionarci,
   !   ! forse serve un buffer da sommare a q alla fine del ciclo su coils number)
   !   q(VAR_JX,:,:,:,:) = 0._R8P
   !   q(VAR_JY,:,:,:,:) = 0._R8P
   !   q(VAR_JZ,:,:,:,:) = 0._R8P
   !   g = 10._R8P*(time_s/td)**3 - 15._R8P*(time_s/td)**4 + 6._R8P*(time_s/td)**5
   !   do n=1, self%coil%total_coils_number
   !   do b=1, blocks_number
   !   do k=1, nk
   !   do j=1, nj
   !   do i=1, ni
   !      coil_id = n
   !      ! use step function to avoid the following original if
   !      !if (time < td) then
   !      !   current_density = g*A(coil_id)/((d(coil_id)-dx)**2)*cos(phase(coil_id)*pi/180.0_R8P)
   !      !else
   !      !   current_density = A(coil_id)/((d(coil_id)-dx)**2)*cos(2*pi*f(coil_id)*(time-td) + &
   !      !   phase(coil_id)*pi/180.0_R8P)
   !      !endif
   !      w_   = nint(sign(1._R8P,td-time_s) + 1._R8P)/2   ! = 1 if td>time,            = 0                              if td<time
   !      w_c_ = 1_I4P - w_                                ! = 0 if td>time,            = 1                              if td<time
   !      g_   = w_ * g + w_c_                             ! = g if td>time,            = 1                              if td<time
   !      f_   = w_c_ * 2._R8P*PI*f(coil_id)*(time_s-td)   ! = 0 if td>time,            = 2._R8P*PI*f(coil_id)*(time-td) if td<time
   !      current_density = g_ * A(coil_id) * cos(f_ + phase(coil_id)*PI/180.0_R8P)
!
   !      ! Lo tengo qui, ma a pensarci bene dovrebbe andare bene così come abiamo fatto (quella sfasata resta a 0)
   !      !f_   = w_c_ * (2._R8P*PI*f(coil_id)*(time_s-td) + phase(coil_id)*PI/180.0_R8P)
   !                        ! = 0 if td>time, = 2._R8P*PI*f(coil_id)*(time-td)+phase(coil_id)*PI/180.0_R8P if td<time
   !      !current_density = g_ * A(coil_id) * cos(f_)*j_vec(4,i,j,k,b)
!
   !      q(VAR_JX,i,j,k,b) = q(VAR_JX,i,j,k,b) + current_density * J_vec(n,1,i,j,k,b)
   !      q(VAR_JY,i,j,k,b) = q(VAR_JY,i,j,k,b) + current_density * J_vec(n,2,i,j,k,b)
   !      q(VAR_JZ,i,j,k,b) = q(VAR_JZ,i,j,k,b) + current_density * J_vec(n,3,i,j,k,b)
   !   enddo
   !   enddo
   !   enddo
   !   enddo
   !enddo
   !!endif
   !   current_density_o = g_ * A(1) * cos(w_c_*2._R8P*PI*f(1)*(time_s-td)+phase(1)*PI/180.0_R8P)
   !   call write_current_behavior_tab('current_density_coil_1.dat', time=time_s, current_density=current_density_o)
   !   current_density_o = g_ * A(4) * cos(w_c_*2._R8P*PI*f(4)*(time_s-td)+phase(4)*PI/180.0_R8P)
   !   !print *, cos(w_c_*2._R8P*PI*f(4)*(time_s-td)+phase(4)*PI/180.0_R8P)
   !   !print *, w_c_*2._R8P*PI*f(4)*(time_s-td)
   !   !print*, phase(4)*PI/180.0_R8P
   !   call write_current_behavior_tab('current_density_coil_4.dat', time=time_s, current_density=current_density_o)
   !endif
   !endassociate
   !endsubroutine compute_coils_current

   subroutine compute_coils_current(self, q, gamma)
   !< Compute current coils sources (DC/AC with smooth envelope).
   class(prism_cpu_object), intent(inout)        :: self
   real(R8P),               intent(inout)        :: q(1:,          &
                                                     1-self%ngc:, &
                                                     1-self%ngc:, &
                                                     1-self%ngc:, &
                                                     1:)
   real(R8P),               intent(in), optional :: gamma
   character(len=128)                            :: fname
   real(R8P)                                     :: current_density
   real(R8P)                                     :: current_density_o
   real(R8P)                                     :: time_s
   real(R8P)                                     :: s, g
   real(R8P)                                     :: phi_rad, omega, theta
   real(R8P)                                     :: f_abs
   integer(I4P)                                  :: w_ac  ! =1 AC, =0 DC (branchless)
   integer(I4P)                                  :: coil_id
   integer(I4P)                                  :: i,j,k,b,n
   real(R8P),                          parameter :: f_tol = 1.0e-30_R8P

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, blocks_number=>self%blocks_number, &
             time=>self%time%time, dt=>self%time%dt, td=>self%coil%td,                 &
             A=>self%coil%A, f=>self%coil%f, phase=>self%coil%phase,                   &
             J_vec=>self%coil%J_vec,                                                   &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz)

      if (present(gamma)) then
         time_s = time + dt*gamma
      else
         time_s = time
      end if
      if (self%coil%total_coils_number >= 1_I4P) then

         ! Azzero termini sorgenti (NB: col PIC potresti voler accumulare in un buffer)
         q(var_Jx,:,:,:,:) = 0._R8P
         q(var_Jy,:,:,:,:) = 0._R8P
         q(var_Jz,:,:,:,:) = 0._R8P

         ! Envelope C^2: clamp(s) in [0,1], g(0)=0, g(1)=1, g'(0)=g'(1)=0, g''(0)=g''(1)=0
         if (td > 0._R8P) then
            s = time_s / td
         else
            s = 1._R8P
         end if
         s = max(0._R8P, min(1._R8P, s))
         g = 10._R8P*s**3 - 15._R8P*s**4 + 6._R8P*s**5

         do n=1, self%coil%total_coils_number
            coil_id = n

            phi_rad = phase(coil_id) * PI / 180._R8P
            omega   = 2._R8P * PI * f(coil_id)

            ! Se f ~ 0 -> DC (w_ac=0). Se f != 0 -> AC (w_ac=1). Branchless.
            f_abs = abs(f(coil_id))
            w_ac  = nint( (sign(1._R8P, f_abs - f_tol) + 1._R8P) * 0.5_R8P )

            ! Theta: DC -> theta = phi ; AC -> theta = omega*(t-td) + phi
            theta = w_ac * omega * (time_s - td) + phi_rad

            ! Unica formula: DC e AC
            current_density = A(coil_id) * g * cos(theta)

            do b=1, blocks_number
               do k=1, nk
                  do j=1, nj
                     do i=1, ni
                        q(var_Jx,i,j,k,b) = q(var_Jx,i,j,k,b) + current_density * J_vec(1,i,j,k,b,n)
                        q(var_Jy,i,j,k,b) = q(var_Jy,i,j,k,b) + current_density * J_vec(2,i,j,k,b,n)
                        q(var_Jz,i,j,k,b) = q(var_Jz,i,j,k,b) + current_density * J_vec(3,i,j,k,b,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         ! Diagnostiche
         do n = 1, self%coil%total_coils_number
            theta = nint( (sign(1._R8P, abs(f(n)) - f_tol) + 1._R8P) * 0.5_R8P ) * (2._R8P*PI*f(n)) * (time_s - td) + &
                    phase(n)*PI/180._R8P
            current_density_o = A(n) * g * cos(theta)
            write(fname,'(A,SS,I0,A)') 'current_density_coil_', n, '.dat'
            call write_current_behavior_tab(trim(fname), time=time_s, current_density=current_density_o)
            !call write_current_behavior_tab('current_density_coil_'//trim(adjustl(str(n))), &
            !                                 time=time_s, current_density=current_density_o)
         enddo
      endif
   endassociate
   endsubroutine compute_coils_current

   subroutine apply_fWL_correction(self, q)
   !< Apply correction if a fWL is present
   class(prism_cpu_object), intent(inout) :: self            !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)           !< Conservative variables.
   integer(I4P)                           :: face            !< Counter
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number,         &
            f=>self%fWLayer%f, layer=>self%fWLayer%layer, C=>self%fWLayer%C, ni_fWL=>self%fWLayer%ni_fWL,     &
            nj_fWL=>self%fWLayer%nj_fWL, nk_fWL=>self%fWLayer%nk_fWL, n=>self%fWLayer%n, s2=>self%fWLayer%s2, &
            alfa_D=>self%fWLayer%alfa_D, alfa_B=>self%fWLayer%alfa_B, beta_D=>self%fWLayer%beta_D,            &
            beta_B=>self%fWLayer%beta_B)
   if (C>0) then
      do face=1, 6
         if (layer(face)) call apply_fWL_correction_fun(blocks_number = blocks_number,      &
                                                        ngc           = ngc,                &
                                                        ni1           = ni_fWL(1,face),     &
                                                        ni2           = ni_fWL(2,face),     &
                                                        nj1           = nj_fWL(1,face),     &
                                                        nj2           = nj_fWL(2,face),     &
                                                        nk1           = nk_fWL(1,face),     &
                                                        nk2           = nk_fWL(2,face),     &
                                                        n             = n(face),            &
                                                        s2            = s2(face),           &
                                                        alfa_D        = alfa_D(face),       &
                                                        beta_D        = beta_D(face),       &
                                                        alfa_B        = alfa_B(face),       &
                                                        beta_B        = beta_B(face),       &
                                                        f             = self%fWLayer%f,     &
                                                        q             = q)
      enddo
   endif
   endassociate
   endsubroutine apply_fWL_correction

   subroutine set_boundary_conditions(self, q, s)
   !< Set boundary conditions of equation.
   class(prism_cpu_object), intent(inout) :: self                 !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:)    !< Conservative variables.
   integer(I4P),  optional, intent(in)    :: s !< Stage counter.
   integer(I4P)                           :: b, c, i, j, k, v        !< Counter.
   integer(I4P)                           :: idelta,jdelta,kdelta    !< IJK delta step for extrapolation.
   integer(I4P)                           :: idelta_n, jdelta_n, kdelta_n !< IJK delta step for Neumann BC.
   integer(I4P)                           :: bc_type                 !< Boundary condition type.
   integer(I4P)                           :: crown                   !< Crown counter.
   integer(I4P)                           :: fec                     !< Boundary fec (1 to 26).
   integer(I4P)                           :: fec_1_6                 !< Boundary fec (1 to 6).
   integer(I4P)                           :: alfa_D, beta_D, gamma_D !< Indici alfa beta gamma come in Barbas.
   integer(I4P)                           :: alfa_B, beta_B, gamma_B !< Indici alfa beta gamma come in Barbas.
   real(R8P)                              :: s1                      !< Coefficiente pari a +-1.
   real(R8P)                              :: ds                      !< Distanza tra le celle in x, y o z.
   real(R8P)                              :: ngc_r, crown_r          !< Numero di gc totale, reale
   real(R8P)                              :: ref(1:9)                !< Vettore di stato di riferimento per assegnazione gc.

   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown,                                                             &
             nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:),               &
             dz=>self%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, dt=>self%time%dt, chi=>self%physics%chi,           &
             nv_c=>self%physics%nv_c, nv_cl=>self%physics%nv_cl, constrained_transport_B=>self%numerics%constrained_transport_B, &
             constrained_transport_D=>self%numerics%constrained_transport_D)

   if (allocated(self%field%maps%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(local_map_bc_crown, dim=1)
            b = local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc_crown(c, 2 ,crown)
               j       = local_map_bc_crown(c, 3 ,crown)
               k       = local_map_bc_crown(c, 4 ,crown)
               idelta  = local_map_bc_crown(c, 5 ,crown)
               jdelta  = local_map_bc_crown(c, 6 ,crown)
               kdelta  = local_map_bc_crown(c, 7 ,crown)
               bc_type = local_map_bc_crown(c, 8 ,crown)
               fec     = local_map_bc_crown(c, 9 ,crown) !da qua la faccia e quindi la normale
               fec_1_6 = fec_1_6_array(fec)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv!(nv_c-nv_cl)
                     q(v,i,j,k,b) = q(v,i-idelta,j-jdelta,k-kdelta,b) !ni,j,k coordinate della cella da cui prendo i valori
                  enddo
               elseif (bc_type == BC_NEUMANN) then
                  if (fec == 1) then
                     idelta_n = nint(abs(real(i))*idelta,kind=I4P) - 1_I4P
                     jdelta_n = jdelta
                     kdelta_n = kdelta
                     do v=1, nv!(nv_c-nv_cl)
                        q(v,i,j,k,b) = q(v,-idelta_n,j-jdelta_n,k-kdelta_n,b)
                     enddo
                  elseif (fec == 2) then
                     idelta_n = 2_I4P*nint(abs(real(i-ni))*idelta,kind=I4P) - 1_I4P
                     jdelta_n = jdelta
                     kdelta_n = kdelta
                     do v=1, nv!(nv_c-nv_cl)
                        q(v,i,j,k,b) = q(v,i-idelta_n,j-jdelta_n,k-kdelta_n,b)
                     enddo
                  elseif (fec == 3) then
                     idelta_n = idelta
                     jdelta_n = nint(abs(real(j))*jdelta,kind=I4P) - 1_I4P
                     kdelta_n = kdelta
                     do v=1, nv!(nv_c-nv_cl)
                        q(v,i,j,k,b) = q(v,i-idelta_n,-jdelta_n,k-kdelta_n,b)
                     enddo
                  elseif (fec == 4) then
                     idelta_n = idelta
                     jdelta_n = 2_I4P*nint(abs(real(j-nj))*jdelta,kind=I4P) - 1_I4P
                     kdelta_n = kdelta
                     do v=1, nv!(nv_c-nv_cl)
                        q(v,i,j,k,b) = q(v,i-idelta_n,j-jdelta_n,k-kdelta_n,b)
                     enddo
                  elseif (fec == 5) then
                     idelta_n = idelta
                     jdelta_n = jdelta
                     kdelta_n = nint(abs(real(k))*kdelta,kind=I4P) - 1_I4P
                     do v=1, nv!(nv_c-nv_cl)
                        q(v,i,j,k,b) = q(v,i-idelta_n,j-jdelta_n,-kdelta_n,b)
                     enddo
                  elseif (fec == 6) then
                     idelta_n = idelta
                     jdelta_n = jdelta
                     kdelta_n = 2_I4P*nint(abs(real(k-nk))*kdelta,kind=I4P) - 1_I4P
                     do v=1, nv!(nv_c-nv_cl)
                        q(v,i,j,k,b) = q(v,i-idelta_n,j-jdelta_n,k-kdelta_n,b)
                     enddo
                  endif
               elseif (bc_type == BC_Silver_Muller) then !Al momento scritta per funzionare solo con un secondo ordine
                  !print *, fec
                  if (fec <= 6) then
                     select case(fec)
                     !Identifico gli alfa beta gamma come nel paper di Barbas, distinguendo tra alfa_D e alfa_B ecc

                     case(1) ! x- face alfa = 2, beta = 3, gamma = 1
                        s1 = -1.0_R8P
                        alfa_D = 2_I4P
                        beta_D = 3_I4P
                        gamma_D = 1_I4P
                        alfa_B = 5_I4P
                        beta_B = 6_I4P
                        gamma_B = 4_I4P
                        ds = dx(b) !distanza tra le celle in x
                        ref = q(:,1,j,k,b) !vettore di stato di riferimento per assegnazione gc

                     case(2) ! x+ face
                        s1 = 1.0_R8P
                        alfa_D = 2_I4P
                        beta_D = 3_I4P
                        gamma_D = 1_I4P
                        alfa_B = 5_I4P
                        beta_B = 6_I4P
                        gamma_B = 4_I4P
                        ref = q(:,ni,j,k,b)

                     case(3) ! y- face
                        s1 = -1.0_R8P
                        alfa_D = 3_I4P
                        beta_D = 1_I4P
                        gamma_D = 2_I4P
                        alfa_B = 6_I4P
                        beta_B = 4_I4P
                        gamma_B = 5_I4P
                        ref = q(:,i,1,k,b)

                     case(4) ! y+ face
                        s1 = 1.0_R8P
                        alfa_D = 3_I4P
                        beta_D = 1_I4P
                        gamma_D = 2_I4P
                        alfa_B = 6_I4P
                        beta_B = 4_I4P
                        gamma_B = 5_I4P
                        ref = q(:,i,nj,k,b)

                     case(5) ! z- face
                        s1 = -1.0_R8P
                        alfa_D = 1_I4P
                        beta_D = 2_I4P
                        gamma_D = 3_I4P
                        alfa_B = 4_I4P
                        beta_B = 5_I4P
                        gamma_B = 6_I4P
                        ref = q(:,i,j,1,b)

                     case(6) ! z+ face
                        s1 = 1.0_R8P
                        alfa_D = 1_I4P
                        beta_D = 2_I4P
                        gamma_D = 3_I4P
                        alfa_B = 4_I4P
                        beta_B = 5_I4P
                        gamma_B = 6_I4P
                        ref = q(:,i,j,nk,b)
                     endselect
                     q(alfa_D,i,j,k,b)  = s1*C0*ref(beta_B)*EPS0
                     q(beta_D,i,j,k,b)  = -s1*C0*ref(alfa_B)*EPS0
                     q(gamma_D,i,j,k,b) = ref(gamma_D)
                     q(alfa_B,i,j,k,b)  = -s1/C0*ref(beta_D)/EPS0
                     q(beta_B,i,j,k,b)  = s1/C0*ref(alfa_D)/EPS0
                     q(gamma_B,i,j,k,b) = ref(gamma_B)

                     do v=(nv_c-nv_cl+1), nv
                        q(v,i,j,k,b) = q(v,i-idelta,j-jdelta,k-kdelta,b)
                     enddo

                  endif
               elseif (bc_type == BC_DIRICHLET) then
                     do v=1, nv
                        q(v,i,j,k,b) = 0._R8P
                     enddo
               elseif (bc_type == BC_PERIOD) then
                  q(:,i,j,k,b) = 0._R8P
                  select case(fec_1_6)
                  case(1)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,ni+i,j   ,k   ,b)
                  case(2)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i-ni,j   ,k   ,b)
                  case(3)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,nj+j,k   ,b)
                  case(4)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,j-nj,k   ,b)
                  case(5)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,j   ,nk+k,b)
                  case(6)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,j   ,k-nk,b)
                  endselect
               endif
            endif
         enddo
      enddo
   endif

   if (self%bc%bc_type(1) == BC_radiative .or. self%bc%bc_type(2) == BC_radiative &
       .or. self%bc%bc_type(3) == BC_radiative .or. self%bc%bc_type(4) == BC_radiative &
       .or. self%bc%bc_type(5) == BC_radiative .or. self%bc%bc_type(6) == BC_radiative) then !Al momento scritta per funzionare solo con un secondo ordine
      if (present(s)) then
         if (s==1_I4P) call self%rk_bc%initialize_stages(q=q)
         if (self%ib%solids_number>0) then !calcolo stadio per le BC
            call self%rk_bc%compute_stage(s=s, dt=self%time%dt, phi=self%ib%phi)
         else
            call self%rk_bc%compute_stage(s=s, dt=self%time%dt)
         endif
         !Calcolo i residui per l'integrazione temporale delle BC (in un futuro da allineare con operatore spaziale qualsiasi)
         call self%compute_residuals_BC(s=s)
         !Imponi effettivamente la BC su q: unico punto del ciclo in cui si "uniscono"
         !Quindi basta cambiare gli indici di quel do per imporlo su una sola variabile, eventualmente
         !(O cambiare i cicli da 1:nv_c a nv_c-nv_cl+1:nv_c)
         if (allocated(self%field%maps%local_map_bc_crown)) then
            do crown=1, ngc
               do c=1, size(local_map_bc_crown, dim=1)
                  b = local_map_bc_crown(c, 1 ,crown)
                  if (b>0) then
                     i       = local_map_bc_crown(c, 2 ,crown)
                     j       = local_map_bc_crown(c, 3 ,crown)
                     k       = local_map_bc_crown(c, 4 ,crown)
                     idelta  = local_map_bc_crown(c, 5 ,crown)
                     jdelta  = local_map_bc_crown(c, 6 ,crown)
                     kdelta  = local_map_bc_crown(c, 7 ,crown)
                     bc_type = local_map_bc_crown(c, 8 ,crown)
                     fec     = local_map_bc_crown(c, 9 ,crown) !da qua la faccia e quindi la normale
                     fec_1_6 = fec_1_6_array(fec)
                     if (bc_type == BC_radiative) then
                        do v=1, nv_c
                           q(v,i,j,k,b) = 2*self%rk_bc%q_bc_rk(v,i,j,k,b,s)-q(v,i-idelta,j-jdelta,k-kdelta,b)
                        enddo
                     endif
                  endif
               enddo
            enddo
         endif
         !Concludi assegnando lo stadio
         if (self%ib%solids_number>0) then
            call self%rk_bc%assign_stage(s=s, phi=self%ib%phi)
         else
            call self%rk_bc%assign_stage(s=s)
         endif
      else !Mi serve solo per il t0, tanto ic è il vuoto praticamente sempre
         q(v,i,j,k,b) = 0.0_R8P
      endif
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine compute_residuals_BC(self,s)
   !< Compute residuals BCs.
   !< La sua scrittura si lega all'ordine di interpolazione dell'operatore spaziale. Al momento
   !< e' scritto per operatore di secondo ordine (1 punto ghost).
   class(prism_cpu_object), intent(inout) :: self                    !< The equation.
   integer(I4P),            intent(in)    :: s                       !< Stage counter.
   real(R8P)                              :: ds                      !< Distanza tra le celle in x, y o z.
   integer(I4P)                           :: i,j,k,b,c,v             !< Counter
   integer(I4P)                           :: idelta,jdelta,kdelta    !< IJK delta step for extrapolation.
   integer(I4P)                           :: bc_type                 !< Boundary condition type.
   integer(I4P)                           :: crown                   !< Crown counter.
   integer(I4P)                           :: fec                     !< Boundary fec (1 to 26).
   integer(I4P)                           :: fec_1_6                 !< Boundary fec (1 to 6).
   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown,                                                     &
             nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:),       &
             dz=>self%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, dt=>self%time%dt, chi=>self%physics%chi,   &
             nv_c=>self%physics%nv_c, nv_cl=>self%physics%nv_cl, div_corr_var=>self%numerics%div_corr_var,               &
             constrained_transport_B=>self%numerics%constrained_transport_B,                                             &
             constrained_transport_D=>self%numerics%constrained_transport_D, q_rk=>self%rk%q_rk,                         &
             q_bc_rk=>self%rk_bc%q_bc_rk,dq_bc_rk=>self%rk_bc%dq_bc_rk)
   if (allocated(self%field%maps%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(local_map_bc_crown, dim=1)
            b = local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               bc_type = local_map_bc_crown(c, 8 ,crown)
               i       = local_map_bc_crown(c, 2 ,crown)
               j       = local_map_bc_crown(c, 3 ,crown)
               k       = local_map_bc_crown(c, 4 ,crown)
               idelta  = local_map_bc_crown(c, 5 ,crown)
               jdelta  = local_map_bc_crown(c, 6 ,crown)
               kdelta  = local_map_bc_crown(c, 7 ,crown)
               fec     = local_map_bc_crown(c, 9 ,crown) !da qua la faccia e quindi la normale
               fec_1_6 = fec_1_6_array(fec)
               if (fec <= 6) then
                  select case(fec)
                  case(1) !xmin
                     ds = dx(b)
                  case(2) !xmax
                     ds = dx(b)
                  case(3) !ymin
                     ds = dy(b)
                  case(4) !ymax
                     ds = dy(b)
                  case(5) !zmin
                     ds = dz(b)
                  case(6) !zmax
                     ds = dz(b)
                  end select
                  do v = 1,6
                     dq_bc_rk(v,i,j,k,b) = -C0*(q_bc_rk(v,i,j,k,b,s)-q_rk(v,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  enddo
                  if (nv_cl == 1_I4P) then
                     dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
                                                q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  elseif (nv_cl == 2_I4P) then
                     dq_bc_rk(nv_c-1,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c-1,i,j,k,b,s)- &
                                                q_rk(nv_c-1,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                     dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
                                                q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  endif
                  !if (div_corr_var == DIV_CORR_VAR_HYPER) then
                  !   if (constrained_transport_D .and. .not.constrained_transport_B) &
                  !      dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
                  !                                 q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  !   if (.not.constrained_transport_D .and. constrained_transport_B) &
                  !      dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
                  !                                 q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  !   if (constrained_transport_D .and. constrained_transport_B) &
                  !      dq_bc_rk(nv_c-1,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c-1,i,j,k,b,s)- &
                  !                                 q_rk(nv_c-1,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  !   if (constrained_transport_D .and. constrained_transport_B) &
                  !      dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
                  !                                 q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
                  !endif
               endif
            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine compute_residuals_BC

   subroutine set_initial_conditions(self) !DA CORREGGERE CON NV_PIC QUANDO SERVE PER BC CARICA SE MODELLO PIC ATTIVO
   !< Set initial conditions and coils on field.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field, q=self%q)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call self%particle_injection%set_particle_initial_injection(field=self%field, pic=self%pic, q_pic=self%q_pic)
      call write_initial_injection_tab(filename='particle_injection.dat', q_pic=self%q_pic, np=self%pic%particle_number)
      call write_initial_injection_tab(filename='neighbour_list.dat', q_pic=real(self%pic%neighbour_list,R8P), &
                                       np=self%pic%particle_number)
   endif
   !call self%coil%set_coils(physics=self%physics, field=self%field) !Lo metto dopo perchè l'interpolatore di correnti azzera
                                                                    !tutto per poter poi fare la sommatoria al relativo tempo

   call self%initialize_coils

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call self%pic%current_weighting(field=self%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
      call self%pic%particle_weighting(field=self%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
      call self%pic%field_weighting(field=self%field, q=self%q, q_pic=self%q_pic, pic_fields=self%pic_fields, nv=self%nv)
   endif
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step, s)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(prism_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),               intent(inout)        :: q(1:,         &
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1:)           !< Conservative variables.
   integer(I4P),            intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
   integer(I4P),            intent(in), optional :: s               !< Stage counter.
   logical                                       :: do_local_update !< Flag for triggering local update.
   logical                                       :: do_set_bc       !< Flag for triggering setting bc.

   ! perform local update if step is not speficied or if first step is selected
   do_local_update = .false.
   do_set_bc       = .false.
   if (.not.present(step)) then
      do_local_update = .true.
      do_set_bc       = .true.
   else
      if (step==1) do_local_update = .true.
      if (step==3) do_set_bc       = .true.
   endif
   if (do_local_update) call self%field%update_ghost_local(q=q)
                        call self%field%update_ghost_mpi(q=q, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q=q, s=s)
   endsubroutine update_ghost

   ! numerical methods
   subroutine compute_auxiliary_fields(self)
   !< Compute auxiliary fields.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   if (self%io%save_divergence_fields) then
      call self%compute_divergence(ivar=VAR_DX,q=self%q,divergence=self%divergence(1,:,:,:,:))
      call self%compute_divergence(ivar=VAR_BX,q=self%q,divergence=self%divergence(2,:,:,:,:))
      call self%compute_divergence(ivar=self%physics%var_Jx,q=self%q,divergence=self%divergence(3,:,:,:,:))
   endif
   if (self%io%save_curl_fields) then
      call self%compute_curl(ivar=VAR_DX,q=self%q,curl=self%curl(1:3,:,:,:,:))
      call self%compute_curl(ivar=VAR_BX,q=self%q,curl=self%curl(4:6,:,:,:,:))
      call self%compute_curl(ivar=self%physics%var_Jx,q=self%q,curl=self%curl(7:9,:,:,:,:))
   endif
   endsubroutine compute_auxiliary_fields

   subroutine compute_dt(self)
   class(prism_cpu_object), intent(inout) :: self                            !< The equation.
   real(R8P)                              :: umax                            !< Maximum speed of waves propagation (light speed).
   real(R8P)                              :: dxyz_min                        !< Minimal space step.
   real(R8P)                              :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                           :: b                               !< Counter.

   dxyz_min = huge(1._R8P)
   associate(blocks_number=>self%blocks_number, dxyz=>self%field%dxyz,chi=>self%physics%chi, evmax=>self%physics%evmax)
   call compute_dxyz_min(blocks_number=blocks_number, dxyz=dxyz, dxyz_min=dxyz_min)
   umax = evmax
   self%time%dt = self%time%CFL*dxyz_min / umax
   endassociate
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endsubroutine compute_dt

   subroutine compute_energy(self)
   !< Compute energy.
   class(prism_cpu_object), intent(inout) :: self          !< The equation.
   real(R8P)                              :: energy_D      !< Energy of D field.
   real(R8P)                              :: energy_B      !< Energy of B field.
   real(R8P)                              :: coil_power    !< Coil power.
   real(R8P)                              :: Poynting_flux !< Total Poynting flux from boundary.

   energy_D = 0.0_R8P
   energy_B = 0.0_R8P
   coil_power = 0.0_R8P
   Poynting_flux = 0.0_R8P
   call compute_e(ivar=VAR_DX, energy=energy_D)
   call compute_e(ivar=VAR_BX, energy=energy_B)
   if (self%coil%total_coils_number > 0_I4P) then
      call compute_coil_power(ivar=self%physics%var_Jx, coil_power=coil_power)
   endif
   call compute_Poynting_flux(Poynting_flux=Poynting_flux)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_B, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, coil_power, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, Poynting_flux, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
   if (allocated(self%energy_D).and.allocated(self%energy_B) &
      .and.allocated(self%coil_power).and.allocated(self%Poynting_flux)) then
      self%energy_D = [self%energy_D, energy_D]
      self%energy_B = [self%energy_B, energy_B]
      self%coil_power = [self%coil_power, coil_power]
      self%Poynting_flux = [self%Poynting_flux, Poynting_flux]
   else
      allocate(self%energy_D(1:self%time%it))
      allocate(self%energy_B(1:self%time%it))
      allocate(self%coil_power(1:self%time%it))
      allocate(self%Poynting_flux(1:self%time%it))
      self%energy_D = energy_D
      self%energy_B = energy_B
      self%coil_power = coil_power
      self%Poynting_flux = Poynting_flux
   endif
   contains
      subroutine compute_e(ivar, energy)
      !< Compute electric/magnetic energy of vector field starting from ivar.
      integer(I4P), intent(in)  :: ivar    !< Starting position of vector field.
      real(R8P),    intent(out) :: energy  !< Energy of the vector field starting from ivar.
      integer(I4P)              :: i,j,k,b !< Counter.
      real(R8P)                 :: cost    !< Costant for the energy computation.

      if (ivar==VAR_DX) then
         cost = EPS0
      elseif (ivar==VAR_BX) then
         cost = MU0
      endif
      energy = 0.0_R8P
      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number, &
                dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:))
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         energy = energy + 0.5_R8P * (self%q(ivar  ,i,j,k,b)*self%q(ivar  ,i,j,k,b) + &
                                      self%q(ivar+1,i,j,k,b)*self%q(ivar+1,i,j,k,b) + &
                                      self%q(ivar+2,i,j,k,b)*self%q(ivar+2,i,j,k,b))/cost*(dx(b)*dy(b)*dz(b))
      enddo
      enddo
      enddo
      enddo
      endassociate
      endsubroutine compute_e

      subroutine compute_coil_power(ivar, coil_power) !A voler essere precisi non è un'energia ma una potenza (per unità di volume)
      !< Compute coil power of vector field starting from ivar.
      integer(I4P), intent(in)  :: ivar        !< Starting position of vector field.
      real(R8P),    intent(out) :: coil_power  !< Coil power of the vector field.
      real(R8P)                 :: cost        !< Costant for the coil power computation.
      integer(I4P)              :: i,j,k,b     !< Counter.

      coil_power = 0.0_R8P
      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number, &
                coil_flag=>self%coil%coil_flag, dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:))
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         if (coil_flag(i,j,k,b) /= 0_I4P) then
            coil_power = coil_power - (self%q(VAR_DX  ,i,j,k,b)*self%q(ivar  ,i,j,k,b) + &
                                       self%q(VAR_DX+1,i,j,k,b)*self%q(ivar+1,i,j,k,b) + &
                                       self%q(VAR_DX+2,i,j,k,b)*self%q(ivar+2,i,j,k,b))/EPS0*(dx(b)*dy(b)*dz(b))
         endif
      enddo
      enddo
      enddo
      enddo
      endassociate
      endsubroutine compute_coil_power

      subroutine compute_Poynting_flux(Poynting_flux)
      !< Compute Poynting flux.
      real(R8P),    intent(out) :: Poynting_flux  !< Power irradiated outside computational domain.
      integer(I4P)              :: i,j,k,b,v      !< Counter.
      real(R8P)                 :: q_boundary(6)  !< Variables at boundary for the Poynting flux computation.
      real(R8P)                 :: n(3)           !< Boundary normal direction

      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number, &
                s=>self%numerics%fdv_half_stencils(1), dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:))
      Poynting_flux = 0.0_R8P
      !Faccia -x
      n = [-1.0_R8P, 0.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,1-s:s,j,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3), q_boundary(4:6))/MU0,n)*(dy(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia +x
      n = [1.0_R8P, 0.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,ni+1-s:ni+s,j,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dy(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia -y
      n = [0.0_R8P, -1.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,1-s:s,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia +y
      n = [0.0_R8P, 1.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,nj+1-s:nj+s,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia -z
      n = [0.0_R8P, 0.0_R8P, -1.0_R8P]
      do b=1, blocks_number
      do j=1, nj
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,j,1-s:s,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dy(b))
      enddo
      enddo
      enddo
      !Faccia +z
      n = [0.0_R8P, 0.0_R8P, 1.0_R8P]
      do b=1, blocks_number
      do j=1, nj
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,j,nk+1-s:nk+s,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dy(b))
      enddo
      enddo
      enddo
      endassociate
   endsubroutine compute_Poynting_flux
   endsubroutine compute_energy

   subroutine compute_energy_error(self)
   !< Compute energy error.
   class(prism_cpu_object), intent(inout) :: self       !< The equation.
   real(R8P)                              :: energy_D0  !< Initial energy of D field.
   real(R8P)                              :: energy_B0  !< Initial energy of B field.
   real(R8P), allocatable                 :: error_D(:) !< Error (normalized) of energy of D field.
   real(R8P), allocatable                 :: error_B(:) !< Error (normalized) of energy of B field.

   if (allocated(self%energy_D).and.allocated(self%energy_B)) then
      energy_D0 = self%energy_D(1)
      energy_B0 = self%energy_B(1)

      error_D = abs(self%energy_D - energy_D0) / abs(energy_D0)
      error_B = abs(self%energy_B - energy_B0) / abs(energy_B0)

      self%rms_energy_error_D = sqrt(sum(error_D)/size(error_D))
      self%rms_energy_error_B = sqrt(sum(error_B)/size(error_B))
   endif
   endsubroutine compute_energy_error

   subroutine impose_div_free(self)
   !< Impose divergence-free property.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   associate(constrained_transport_D=>self%numerics%constrained_transport_D,&
             constrained_transport_B=>self%numerics%constrained_transport_B,div_corr_var=>self%numerics%div_corr_var)
   if (constrained_transport_D.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=1_I4P)
   if (constrained_transport_B.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=4_I4P)
   ! here should go also other corrections...
   endassociate
   endsubroutine impose_div_free

   subroutine impose_ct_correction(self, ivar)
   !< Impose Constrained Transport Correction on vectorial variable q(ivar:ivar+2).
   !< Note that self%divergence memory is used as buffer, be careful.
   class(prism_cpu_object), intent(inout) :: self      !< The equation.
   integer(I4P),            intent(in)    :: ivar      !< Variable (start) index in q.
   real(R8P)                              :: dq_max    !< Maximum residual.
   integer(I4P)                           :: iter      !< Counter.
   integer(I4P)                           :: i,j,k,b,v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number, buffer=>self%divergence)
   call self%compute_divergence(ivar=ivar,q=self%q,divergence=buffer(4,:,:,:,:))
   if (blocks_number>0) then
      do iter=1, self%flail%iterations
         call compute_smoothing_multigrid(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=1_I4P,blocks_number=blocks_number, &
                                          dxyz=self%field%dxyz,                                           &
                                          f=-buffer(4:4,:,:,:,:),                                         &
                                          q=buffer(7:7,:,:,:,:),                                          &
                                          dq_max=dq_max,                                                  &
                                          dq=buffer(5:5,:,:,:,:),                                         &
                                          iterations_init=self%flail%iterations_init,                     &
                                          iterations_fine=self%flail%iterations_fine,                     &
                                          iterations_coarse=self%flail%iterations_coarse)
         if (dq_max < self%flail%tolerance) exit
      enddo
      call self%mpih%print_message('FLAIL convergence reached at iteration '//trim(str(iter,.true.)))
      call self%compute_gradient(ivar=1,q=buffer(7:7,:,:,:,:),gradient=buffer(4:6,:,:,:,:))
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, 3
                     self%q(ivar+v-1,i,j,k,b) = self%q(ivar+v-1,i,j,k,b) + buffer(3+v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endassociate
   endsubroutine impose_ct_correction

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(prism_cpu_object), intent(inout) :: self             !< The equation.
   character(*),            intent(in)    :: filename         !< Input file name.
   real(R8P)                              :: timing(1:2)      !< Tic toc timing.
   real(R8P)                              :: timing_step(1:2) !< Tic toc timing.
   real(R8P)                              :: F_l(3)           !< Lorentz force for leapfrog preliminary integration
   integer(I4P)                           :: i                !< Counter.

   ! initialization
   call self%initialize(filename=filename)
   if (self%io%restart) then
      call self%mpih%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call self%mpih%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
   else
      call self%mpih%print_message('impose initial conditions start')
      do i=1, self%ic%amr_iterations
         call self%mpih%print_message('  AMR/set IC iteration:'//trim(str(i,.true.)))
         call self%set_initial_conditions
         !if (self%ib%solids_number > 0) call self%compute_phi()
         !call self%amr_update
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
      call self%mpih%print_message('impose initial conditions finish')
   endif
   !if (self%ib%solids_number > 0) call self%compute_phi()
   ! call self%amr_update
   call self%compute_divergence(ivar=1,q=self%q,divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(ivar=4,q=self%q,divergence=self%divergence(2,:,:,:,:))
   call self%compute_divergence(ivar=7,q=self%q,divergence=self%divergence(3,:,:,:,:))
   call self%save_simulation_data
   call self%update_ghost(q=self%q) ! Aggiunto da FN 
   call self%compute_energy
   !call self%save_energy_error(is_to_open=.true.)
   call self%save_energy_history(is_to_open=.true.)
   call self%save_divergence_history(is_to_open=.true.)
   call self%io%open_file_residuals(nv=self%nv)

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      if(self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) then
         ! first time integration done apart with explicit euler scheme to iniziale leapfrog
         call self%leapfrog_pic%assign_step(s=1, q_pic=self%q_pic)
         call self%compute_dt
         !< Pic residual computation
         !Qua ci va il calcolo dei campi nelle posizioni delle particelle se PIC
         !ma devi metterlo nell'inizializzazione per coerenza
         !< Integration of equations
         self%q_pic(1,:) = self%q_pic(1,:) + self%time%dt * self%q_pic(4,:)
         self%q_pic(2,:) = self%q_pic(2,:) + self%time%dt * self%q_pic(5,:)
         self%q_pic(3,:) = self%q_pic(3,:) + self%time%dt * self%q_pic(6,:)
         do i = 1, self%pic%particle_number
            F_l = crossproduct(self%q_pic(4:6,i), self%pic_fields(4:6,i))
            self%q_pic(4,i) = self%q_pic(4,i)+self%time%dt*self%q_pic(7,i)/self%q_pic(8,i)*(self%pic_fields(1,i)+F_l(1))
            self%q_pic(5,i) = self%q_pic(5,i)+self%time%dt*self%q_pic(7,i)/self%q_pic(8,i)*(self%pic_fields(2,i)+F_l(2))
            self%q_pic(6,i) = self%q_pic(6,i)+self%time%dt*self%q_pic(7,i)/self%q_pic(8,i)*(self%pic_fields(3,i)+F_l(3))
            !self%q_pic(4:6,i) = self%q_pic(4:6,i) + self%time%dt * self%q_pic(8,i) / self%q_pic(7,i) * &
            !                  (self%pic_fields(1:3,i) + crossproduct(self%q_pic(4:6,i), self%pic_fields(4:6,i)))
         enddo
      endif
   endif

   if (self%numerics%scheme_time==NUM_SCHEME_TIME_LEAPFROG) then
      ! first time integration done apart with explicit euler scheme to iniziale leapfrog
      call self%leapfrog%assign_step(s=1, q=self%q)
      call self%compute_dt
      !Qua c'era il calcolo delle correnti delle particelle se PIC
      !Ora è nelle condizioni iniziali per coerenza con if legato a se ho pic o meno
      call self%compute_residuals(q=self%q, dq=self%dq)
      self%q = self%q + self%time%dt * self%dq
   endif

   ! integration
   call self%mpih%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_status(file_name='memory_cpu-'//self%mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (mod(self%time%it,self%amr%frequency)==0) then
         call self%mpih%barrier(tictoc=.true.)
         !call self%amr_update
         call self%mpih%barrier(tictoc=.true.)
      endif

      call self%compute_dt
      if ((self%time%it_max <= 0).and.(self%time%time+self%time%dt > self%time%time_max)) &
         self%time%dt=self%time%time_max-self%time%time

      call self%integrate

      self%time%time = self%time%time + self%time%dt
      call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)

      call self%save_simulation_data
      call self%update_ghost(q=self%q) ! Aggiunto da FN 
      call self%compute_energy
      !call self%save_energy_error
      call self%save_energy_history
      call self%compute_divergence(ivar=1,q=self%q,divergence=self%divergence(1,:,:,:,:))
      call self%compute_divergence(ivar=4,q=self%q,divergence=self%divergence(2,:,:,:,:))
      call self%compute_divergence(ivar=7,q=self%q,divergence=self%divergence(3,:,:,:,:))
      call self%save_divergence_history

      if (((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
         ((self%time%it>=self%time%it_max).and.(self%time%it_max > 0))) exit integration

      call self%mpih%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
   enddo integration
   call self%mpih%barrier(tictoc=.true., timing=timing(2), single=.true.)
   !call self%compute_energy_error
   call self%save_simulation_data
   call self%io%close_file_residuals
   !call self%save_energy_error(is_to_close=.true.)
   !call self%mpih%print_message('Initial/final energy of D field: '//trim(str(sqrt(self%energy_D(1))))//' '//&
   !                                                                  trim(str(sqrt(self%energy_D(size(self%energy_D))))))
   !call self%mpih%print_message('Initial/final energy of B field: '//trim(str(sqrt(self%energy_B(1))))//' '//&
   !                                                                  trim(str(sqrt(self%energy_D(size(self%energy_B))))))
   !call self%mpih%print_message('RMS Error of D field: '//trim(str(self%rms_energy_error_D)))
   !call self%mpih%print_message('RMS Error of B field: '//trim(str(self%rms_energy_error_B)))
   call self%save_energy_history(is_to_close=.true.)
   call self%update_ghost(q=self%q) ! Aggiunto da FN 
   call self%compute_divergence(ivar=1,q=self%q,divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(ivar=4,q=self%q,divergence=self%divergence(2,:,:,:,:))
   call self%compute_divergence(ivar=7,q=self%q,divergence=self%divergence(3,:,:,:,:))
   call self%save_divergence_history(is_to_close=.true.)
   call self%mpih%finalize
   endsubroutine simulate

   ! pointer TBP concrete implementations
   subroutine compute_residuals_fd_centered(self, q, dq, s)
   !< Compute residuals of equation, space operator, centered finite difference schemes.
   class(prism_cpu_object), intent(inout) :: self               !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)              !< Conservative variables.
   real(R8P),               intent(inout) :: dq(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:)             !< Residuals.
   integer(I4P),  optional, intent(in)    :: s !< Stage counter.
   integer(I4P)                           :: i,j,k,b            !< Counter
   real(R8P)                              :: curlD(3), curlB(3) !< Residuals components.
   real(R8P)                              :: KO_Dx_x,KO_Dx_y,KO_Dx_z
   real(R8P)                              :: KO_Dy_x,KO_Dy_y,KO_Dy_z
   real(R8P)                              :: KO_Dz_x,KO_Dz_y,KO_Dz_z
   real(R8P)                              :: KO_Bx_x,KO_Bx_y,KO_Bx_z
   real(R8P)                              :: KO_By_x,KO_By_y,KO_By_z
   real(R8P)                              :: KO_Bz_x,KO_Bz_y,KO_Bz_z
   real(R8P), parameter :: sigma = 1000.01_R8P

   call self%apply_fWL_correction(q=q)
   call self%update_ghost(q=q, s=s)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dxyz=>self%field%dxyz,                                                                                   &
             s1=>self%numerics%fdv_half_stencils(1),                                                                  &
             s4=>self%numerics%fdv_half_stencils(4),                                                                  &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz)
   if (blocks_number > 0) then
      ! compute RHS dD/dt = curl(B/MU0) - J, dB/dt = -curl(D/EPS0)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                       q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),&
                                       curl=curlD)
         call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                       q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),&
                                       curl=curlB)

  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(1,b),q=q(VAR_DX,i-s4:i+s4,j,k,b),d4q_ds4=KO_Dx_x);KO_Dx_x=dxyz(1,b)**3*KO_Dx_x
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(2,b),q=q(VAR_DX,i,j-s4:j+s4,k,b),d4q_ds4=KO_Dx_y);KO_Dx_y=dxyz(2,b)**3*KO_Dx_y
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(3,b),q=q(VAR_DX,i,j,k-s4:k+s4,b),d4q_ds4=KO_Dx_z);KO_Dx_z=dxyz(3,b)**3*KO_Dx_z
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(1,b),q=q(VAR_DY,i-s4:i+s4,j,k,b),d4q_ds4=KO_Dy_x);KO_Dy_x=dxyz(1,b)**3*KO_Dy_x
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(2,b),q=q(VAR_DY,i,j-s4:j+s4,k,b),d4q_ds4=KO_Dy_y);KO_Dy_y=dxyz(2,b)**3*KO_Dy_y
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(3,b),q=q(VAR_DY,i,j,k-s4:k+s4,b),d4q_ds4=KO_Dy_z);KO_Dy_z=dxyz(3,b)**3*KO_Dy_z
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(1,b),q=q(VAR_DZ,i-s4:i+s4,j,k,b),d4q_ds4=KO_Dz_x);KO_Dz_x=dxyz(1,b)**3*KO_Dz_x
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(2,b),q=q(VAR_DZ,i,j-s4:j+s4,k,b),d4q_ds4=KO_Dz_y);KO_Dz_y=dxyz(2,b)**3*KO_Dz_y
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(3,b),q=q(VAR_DZ,i,j,k-s4:k+s4,b),d4q_ds4=KO_Dz_z);KO_Dz_z=dxyz(3,b)**3*KO_Dz_z
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(1,b),q=q(VAR_BX,i-s4:i+s4,j,k,b),d4q_ds4=KO_Bx_x);KO_Bx_x=dxyz(1,b)**3*KO_Bx_x
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(2,b),q=q(VAR_BX,i,j-s4:j+s4,k,b),d4q_ds4=KO_Bx_y);KO_Bx_y=dxyz(2,b)**3*KO_Bx_y
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(3,b),q=q(VAR_BX,i,j,k-s4:k+s4,b),d4q_ds4=KO_Bx_z);KO_Bx_z=dxyz(3,b)**3*KO_Bx_z
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(1,b),q=q(VAR_BY,i-s4:i+s4,j,k,b),d4q_ds4=KO_By_x);KO_By_x=dxyz(1,b)**3*KO_By_x
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(2,b),q=q(VAR_BY,i,j-s4:j+s4,k,b),d4q_ds4=KO_By_y);KO_By_y=dxyz(2,b)**3*KO_By_y
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(3,b),q=q(VAR_BY,i,j,k-s4:k+s4,b),d4q_ds4=KO_By_z);KO_By_z=dxyz(3,b)**3*KO_By_z
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(1,b),q=q(VAR_BZ,i-s4:i+s4,j,k,b),d4q_ds4=KO_Bz_x);KO_Bz_x=dxyz(1,b)**3*KO_Bz_x
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(2,b),q=q(VAR_BZ,i,j-s4:j+s4,k,b),d4q_ds4=KO_Bz_y);KO_Bz_y=dxyz(2,b)**3*KO_Bz_y
  !call compute_derivative4_fd_centered(s=s4,ds=dxyz(3,b),q=q(VAR_BZ,i,j,k-s4:k+s4,b),d4q_ds4=KO_Bz_z);KO_Bz_z=dxyz(3,b)**3*KO_Bz_z

         dq(VAR_DX,i,j,k,b) =  curlB(1)/MU0 - q(var_Jx,i,j,k,b)!- sigma*C0*(KO_Dx_x+KO_Dx_y+KO_Dx_z)/16._R8P
         dq(VAR_DY,i,j,k,b) =  curlB(2)/MU0 - q(var_Jy,i,j,k,b)!- sigma*C0*(KO_Dy_x+KO_Dy_y+KO_Dy_z)/16._R8P
         dq(VAR_DZ,i,j,k,b) =  curlB(3)/MU0 - q(var_Jz,i,j,k,b)!- sigma*C0*(KO_Dz_x+KO_Dz_y+KO_Dz_z)/16._R8P
         dq(VAR_BX,i,j,k,b) = -curlD(1)/EPS0                   !- sigma*C0*(KO_Bx_x+KO_Bx_y+KO_Bx_z)/16._R8P
         dq(VAR_BY,i,j,k,b) = -curlD(2)/EPS0                   !- sigma*C0*(KO_By_x+KO_By_y+KO_By_z)/16._R8P
         dq(VAR_BZ,i,j,k,b) = -curlD(3)/EPS0                   !- sigma*C0*(KO_Bz_x+KO_Bz_y+KO_Bz_z)/16._R8P
      enddo
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine compute_residuals_fd_centered

   subroutine compute_residuals_fv_centered(self, q, dq, s)
   !< Compute residuals of equation, space operator, centered finite volume schemes.
   class(prism_cpu_object), intent(inout) :: self                                             !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)                                            !< Conservative variables.
   real(R8P),               intent(inout) :: dq(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:)                                           !< Residuals.
   integer(I4P),  optional, intent(in)    :: s !< Stage counter.
   integer(I4P)                           :: i,j,k,b,d,v                                      !< Counter
   real(R8P),    parameter                :: sir(3,3) = reshape([1._R8P,0._R8P,0._R8P,&
                                                                 0._R8P,1._R8P,0._R8P,&
                                                                 0._R8P,0._R8P,1._R8P],[3,3]) !< Direction versor, real.

   call self%apply_fWL_correction(q=q)
   call self%update_ghost(q=q)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dxyz=>self%field%dxyz, flxyz_c=>self%flxyz_c, flx_f=>self%flx_f, fly_f=>self%fly_f, flz_f=>self%flz_f,   &
             s=>self%numerics%fdv_half_stencils(1),                                                                   &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz, chi=>self%physics%chi)
   if (blocks_number > 0) then
      ! compute fluxes at cell centers
      do b=1, blocks_number
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
         do d=1,3 ! x, y, z
            call compute_fluxes_Maxwell(sir=sir(:,d),q=q(:,i,j,k,b),f=flxyz_c(:,1,d,i,j,k,b),chi=chi)
         enddo
      enddo
      enddo
      enddo
      enddo
      ! reconstruct fluxes at cell faces
      d = 1
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=0, ni
         do v=1, nv_c
            call compute_reconstruction_r_fv_centered(s=s,q=flxyz_c(v,1,d,1+i-s:i+s,j,k,b),qr=flx_f(v,i,j,k,b))
         enddo
      enddo
      enddo
      enddo
      enddo
      d = 2
      do b=1, blocks_number
      do k=1, nk
      do j=0, nj
      do i=1, ni
         do v=1, nv_c
            call compute_reconstruction_r_fv_centered(s=s,q=flxyz_c(v,1,d,i,1+j-s:j+s,k,b),qr=fly_f(v,i,j,k,b))
         enddo
      enddo
      enddo
      enddo
      enddo
      d = 3
      do b=1, blocks_number
      do k=0, nk
      do j=1, nj
      do i=1, ni
         do v=1, nv_c
            call compute_reconstruction_r_fv_centered(s=s,q=flxyz_c(v,1,d,i,j,1+k-s:k+s,b),qr=flz_f(v,i,j,k,b))
         enddo
      enddo
      enddo
      enddo
      enddo
      ! compute fluxes difference for RHS, dF/ds = (F(i+1/2)-F(i-1/2))/Ds
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         do v=1, nv_c
            dq(v,i,j,k,b) = - (flx_f(v,i,j,k,b) - flx_f(v,i-1,j,k,b)) / dxyz(1,b) &
                            - (fly_f(v,i,j,k,b) - fly_f(v,i,j-1,k,b)) / dxyz(2,b) &
                            - (flz_f(v,i,j,k,b) - flz_f(v,i,j,k-1,b)) / dxyz(3,b)
         enddo
         ! J sources
         dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) - q(var_Jx,i,j,k,b)
         dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) - q(var_Jy,i,j,k,b)
         dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) - q(var_Jz,i,j,k,b)
      enddo
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine compute_residuals_fv_centered

   subroutine compute_residuals_weno(self, q, dq, s)
   !< Compute residuals of equation, space operator, WENO schemes.
   class(prism_cpu_object), intent(inout) :: self   !< The equation.
   real(R8P),               intent(inout) :: q(1:,       &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)  !< Conservative variables.
   real(R8P),               intent(inout) :: dq(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:) !< Residuals.
   integer(I4P),  optional, intent(in)    :: s !< Stage counter.

   call self%apply_fWL_correction(q=q)
   call self%update_ghost(q=q)
   !call self%integrate_eikonal_coils(q=q)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, nv_c=>self%nv_c,blocks_number=>self%blocks_number,&
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:),                                        &
             flx=>self%flxyz_c(:,1,1,:,:,:,:), fly=>self%flxyz_c(:,1,2,:,:,:,:), flz=>self%flxyz_c(:,1,3,:,:,:,:),                &
             weno_s=>self%weno%S, weno_zeps=>self%weno%zeps,                                                                      &
             weno_a=>self%weno%a, weno_p=>self%weno%p, weno_d=>self%weno%d, weno_c=>self%weno%c,                                  &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz, chi=>self%physics%chi,        &
             evmax=>self%physics%evmax, erw=>self%physics%erw, elw=>self%physics%elw)

   if (blocks_number > 0) then
      call compute_fluxes_convective_weno(dir=1,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,               &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_c=weno_c,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,chi=chi,                                                      &
                                     q=q,fluxes=flx)
      call compute_fluxes_convective_weno(dir=2,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,               &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_c=weno_c,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,chi=chi,                                                      &
                                     q=q,fluxes=fly)
      call compute_fluxes_convective_weno(dir=3,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,               &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_c=weno_c,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,chi=chi,                                                      &
                                     q=q,fluxes=flz)
      call compute_fluxes_difference(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv_c=nv_c, &
                                     var_Jx=var_Jx, var_Jy=var_Jy, var_Jz=var_Jz, &
                                     dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, dq=dq, q=q)
   endif
   endassociate
   endsubroutine compute_residuals_weno

   subroutine integrate_blanesmoan(self)
   !< Integrate equation, time operator, Yoshida RK scheme.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   associate(nc=>self%blanesmoan%nc,a=>self%blanesmoan%a,b=>self%blanesmoan%b)
   call self%compute_coils_current(q=self%q)
   do s=1, nc
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s==1) call self%save_residuals
      self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
      self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
      self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
      call self%compute_residuals(q=self%q, dq=self%dq)
      self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
      self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
      self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   enddo
   call self%impose_div_free
   endassociate
   endsubroutine integrate_blanesmoan

   subroutine integrate_cfm(self)
   !< Integrate equation, time operator, Commutator-Free Magnus integrator.
   class(prism_cpu_object), intent(inout) :: self             !< The equation.
   real(R8P), parameter                   :: toll=1.0e-14_R8P !< CFM coefficients tollerance.
   integer(I4P)                           :: s,ss             !< Counter.

   call self%compute_coils_current(q=self%q)
   associate(dt=>self%time%dt,s_coeffs=>self%cfm%s_coeffs,e_coeffs=>self%cfm%e_coeffs)
   self%cfm%q = self%q
   call self%compute_residuals(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,1))
   do s=2, self%cfm%n_stages
      do ss=1, s-1
         if (abs(s_coeffs(s,ss))>toll) &
            call self%cfm%compute_exponential_update(alpha=dt*s_coeffs(s,ss),dq=self%cfm%dq(:,:,:,:,:,ss),q=self%cfm%q)
      enddo
      call self%compute_residuals(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,s))
   enddo
   self%cfm%q = self%q
   do s=1, self%cfm%n_stages
      if (abs(self%cfm%e_coeffs(s))>toll) &
         call self%cfm%compute_exponential_update(alpha=dt*e_coeffs(s),dq=self%cfm%dq(:,:,:,:,:,s),q=self%cfm%q)
   enddo
   self%q = self%cfm%q
   endassociate
   call self%impose_div_free
   endsubroutine integrate_cfm

   subroutine integrate_leapfrog(self)
   !< Integrate equation, time operator, leapfrog scheme.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%compute_coils_current(q=self%q)
   call self%compute_residuals(q=self%q, dq=self%dq)
   call self%save_residuals
   call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)
   call self%impose_div_free
   endsubroutine integrate_leapfrog

   subroutine integrate_leapfrog_pic(self)
   !< Integrate equation, time operator, leapfrog scheme for particle in cell
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   !< Maxwell source terms computation: particles and coils
   call self%pic%particle_cartesian_grid_index(field=self%field, q_pic=self%q_pic)
   call self%pic%current_weighting(field=self%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
   call self%compute_coils_current(q=self%q)
   !< Maxwell residuals computation
   call self%compute_residuals(q=self%q, dq=self%dq)
   call self%save_residuals
   !< Pic residual computation
   call self%pic%field_weighting(field=self%field, q=self%q, q_pic=self%q_pic, pic_fields=self%pic_fields, nv=self%nv)
   !< Integration of equations
   call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)
   call self%leapfrog_pic%integrate(dt=self%time%dt, q_pic=self%q_pic, pic_fields=self%pic_fields)
   call self%impose_div_free
   endsubroutine integrate_leapfrog_pic

   subroutine integrate_rk_ls(self)
   !< Integrate equation, time operator, RK classical low storage schemes.
   !< Low storage RK working on q_rk(:,:,:,:,:,1)/q as stages, update q in place.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%compute_coils_current(q=self%q) !da modificare per avere i tempi corretti
   call self%rk%initialize_stages(q=self%q)
   do s=1, self%rk%nrk
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk%compute_stage_ls(s=s,dt=self%time%dt,phi=self%ib%phi,dq=self%dq,q=self%q)
      else
         call self%rk%compute_stage_ls(s=s,dt=self%time%dt,dq=self%dq,q=self%q)
      endif
   enddo
   call self%impose_div_free
   endsubroutine integrate_rk_ls

   subroutine integrate_rk_ssp(self)
   !< Integrate equation, time operator, SSP RK schemes.
   !< SSP RK working on q_rk as stages.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%sub_external_fields(field=self%field, time=self%time%time, dt=self%time%dt, q=self%q)
   call self%rk%initialize_stages(q=self%q)
   do s=1, self%rk%nrk
      if (self%ib%solids_number>0) then
         call self%rk%compute_stage(s=s, dt=self%time%dt, phi=self%ib%phi)
      else
         call self%rk%compute_stage(s=s, dt=self%time%dt)
      endif
      call self%compute_coils_current(q=self%rk%q_rk(:,:,:,:,:,s), gamma=self%rk%gamm(s))
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq, s=s)
      !if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk%assign_stage(s=s, q=self%dq, phi=self%ib%phi)
      else
         call self%rk%assign_stage(s=s, q=self%dq)
      endif
   enddo
   if (self%ib%solids_number>0) then
      call self%rk%update_q(dt=self%time%dt, phi=self%ib%phi, q=self%q)
      call self%update_q_BC(dt=self%time%dt, phi=self%ib%phi)
   else
      call self%rk%update_q(dt=self%time%dt, q=self%q, dq=self%dq)
      call self%update_q_BC(dt=self%time%dt)
      call self%save_residuals
   endif
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%add_external_fields(field=self%field, time=self%time%time, dt=self%time%dt, q=self%q)
   endsubroutine integrate_rk_ssp

   subroutine integrate_rk_ssp_pic(self)
   !< Integrate equation, time operator, SSP RK schemes.
   !< SSP RK working on q_rk as stages.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   !call sub_external_fields(self = self%external_fields, field = self%field, &
   !                        time = self%time%time, dt = self%time%dt, q = self%q)

   !Inizializzo stadi RK per campi e PIC
   call self%rk%initialize_stages(q=self%q)
   call self%rk_pic%initialize_stages(q_pic=self%q_pic)

   do s=1, self%rk%nrk
      !Calcolo stadio RK per campi e PIC
      if (self%ib%solids_number>0) then
         call self%rk%compute_stage(s=s, dt=self%time%dt, phi=self%ib%phi)
      else
         call self%rk%compute_stage(s=s, dt=self%time%dt)
      endif
      call self%rk_pic%compute_stage(s=s, dt=self%time%dt)
      !Calcolo termini sorgente Maxwell da particelle e bobine
      call self%pic%particle_cartesian_grid_index(field=self%field, q_pic=self%rk_pic%q_pic_rk(:,:,s))
      call self%pic%current_weighting(field=self%field, q=self%rk%q_rk(:,:,:,:,:,s), &
                                       q_pic=self%rk_pic%q_pic_rk(:,:,s), nv=self%nv)
      call self%compute_coils_current(q=self%rk%q_rk(:,:,:,:,:,s), gamma=self%rk%gamm(s))
      !Calcolo residui Maxwell
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq, s=s)
      if (s==1) call self%save_residuals
      !Calcolo residui PIC: calcolati direttamente nell'assegnazione dello stadio RK
      !Interpolo quindi i campi (probabilmente è qui che ti conviene sommare e sottrarre i campi esterni)
      call self%pic%field_weighting(field=self%field, q=self%rk%q_rk(:,:,:,:,:,s), &
                                    q_pic=self%rk_pic%q_pic_rk(:,:,s), pic_fields=self%pic_fields, nv=self%nv)
      !Assegno lo stadio RK per campi e PIC
      if (self%ib%solids_number>0) then
         call self%rk%assign_stage(s=s, q=self%dq, phi=self%ib%phi)
      else
         call self%rk%assign_stage(s=s, q=self%dq)
      endif
      call self%rk_pic%assign_stage(s=s, pic_fields=self%pic_fields)
   enddo
   ! Completo l'integrazione temporale
   if (self%ib%solids_number>0) then
      call self%rk%update_q(dt=self%time%dt, phi=self%ib%phi, q=self%q)
      call self%update_q_BC(dt=self%time%dt, phi=self%ib%phi)
   else
      call self%rk%update_q(dt=self%time%dt, q=self%q)
      call self%update_q_BC(dt=self%time%dt)
   endif
   call self%rk_pic%update_q_pic(dt=self%time%dt, q_pic=self%q_pic)
   !Aggiorno i termini sorgente di Maxwell al tempo in cui andrò a plottare i risultati
   call self%impose_div_free
   call self%pic%particle_cartesian_grid_index(field=self%field, q_pic=self%q_pic)
   call self%pic%current_weighting(field=self%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
   call self%pic%particle_weighting(field=self%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
   call self%compute_coils_current(q=self%q)
   !call add_external_fields(self = self%external_fields, field = self%field, &
   !                        time = self%time%time, dt = self%time%dt, q = self%q)
   endsubroutine integrate_rk_ssp_pic

   subroutine update_q_BC(self, dt, phi)
   !< Update RK q ghost cells.
   class(prism_cpu_object), intent(inout) :: self             !< RK object.
   real(R8P),        intent(in)           :: dt               !< Current time step.
   real(R8P),        intent(in), optional :: phi(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)          !< IB distance.
   integer(I4P)                           :: all_solids       !< Last phi index, all solids summary.
   integer(I4P)                           :: i, j, k, b, v, s, c !< Counter.

   integer(I4P)                           :: idelta,jdelta,kdelta    !< IJK delta step for extrapolation.
   integer(I4P)                           :: bc_type                 !< Boundary condition type.
   integer(I4P)                           :: crown                   !< Crown counter.
   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown,                                                       &
                nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:),      &
                dz=>self%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, dt=>self%time%dt, chi=>self%physics%chi,  &
                nv_c=>self%physics%nv_c, nv_cl=>self%physics%nv_cl,                                                        &
                constrained_transport_B=>self%numerics%constrained_transport_B,                                            &
                constrained_transport_D=>self%numerics%constrained_transport_D, nrk=>self%rk_bc%nrk, &
                q_bc_rk=>self%rk_bc%q_bc_rk, blocks_number=>self%blocks_number, beta=>self%rk_bc%beta)

   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(6) default(firstprivate) shared(phi,self)
      do s=1, nrk
         do b=1, blocks_number
            do k=1-ngc, nk+ngc
               do j=1-ngc, nj+ngc
                  do i=1-ngc, ni+ngc
                     !(O cambiare i cicli da 1:nv_c a nv_c-nv_cl+1:nv_c)
                     do v=1, nv_c
                        if (phi(all_solids,i,j,k,b) < 0._R8P) then
                           q_bc_rk(v,i,j,k,b,nrk+1) = q_bc_rk(v,i,j,k,b,nrk+1) + dt * beta(s) * q_bc_rk(1,i,j,k,b,s)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(6) default(firstprivate) shared(self)
      do s=1, nrk
         do b=1, blocks_number
            do k=1-ngc, nk+ngc
               do j=1-ngc, nj+ngc
                  do i=1-ngc, ni+ngc
                     do v=1, nv_c
                        q_bc_rk(v,i,j,k,b,nrk+1) = q_bc_rk(v,i,j,k,b,nrk+1) + dt * beta(s) * q_bc_rk(1,i,j,k,b,s)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   if (allocated(self%field%maps%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(local_map_bc_crown, dim=1)
            b = local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               bc_type = local_map_bc_crown(c, 8 ,crown)
               if (bc_type == BC_radiative) then
                  i       = local_map_bc_crown(c, 2 ,crown)
                  j       = local_map_bc_crown(c, 3 ,crown)
                  k       = local_map_bc_crown(c, 4 ,crown)
                  idelta  = local_map_bc_crown(c, 5 ,crown)
                  jdelta  = local_map_bc_crown(c, 6 ,crown)
                  kdelta  = local_map_bc_crown(c, 7 ,crown)
                  do v=1, nv_c
                     self%q(nv_c,i,j,k,b) = 2*q_bc_rk(1,i,j,k,b,nrk+1)-self%q(nv_c,i-idelta,j-jdelta,k-kdelta,b)
                  enddo
                  !print *, 'Updating BC SM', b, ' cell (', i, ',', j, ',', k, ')'
               endif
            endif
            !Qua ci aggiungi gli altri if a seconda elle variabili su cui vuoi implementare questa BC
         enddo
      enddo
   endif
   endassociate
   endsubroutine update_q_BC

   subroutine integrate_rk_yoshida(self)
   !< Integrate equation, time operator, Yoshida RK scheme.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%compute_coils_current(q=self%q)
   do s=1, self%rk%nrk - 1
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s==1) call self%save_residuals
      self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
      self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
      self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
      call self%compute_residuals(q=self%q, dq=self%dq)
      self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
      self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
      self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   enddo
   call self%compute_residuals(q=self%q, dq=self%dq)
   self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   call self%impose_div_free
   endsubroutine integrate_rk_yoshida

   ! non TBP
   subroutine compute_dxyz_min(blocks_number, dxyz, dxyz_min)
   !< Compute minimum dxyz space step.
   integer(I4P), intent(in)  :: blocks_number !< Number of blocks.
   real(R8P),    intent(in)  :: dxyz(:,:)     !< XYZ space steps.
   real(R8P),    intent(out) :: dxyz_min      !< Minimum space step.
   integer(I4P)              :: b             !< Counter.

   dxyz_min = huge(0._R8P)
   !$omp parallel do shared(dxyz) reduction(min:dxyz_min)
   do b=1, blocks_number
      dxyz_min = min(dxyz_min, dxyz(1,b), dxyz(2,b), dxyz(3,b))
   enddo
   dxyz_min = dxyz_min * 0.5_R8P
   endsubroutine compute_dxyz_min

   subroutine compute_fluxes_convective_weno(dir,blocks_number,ni,nj,nk,ngc,nv_c,weno_s,weno_a,weno_p,weno_d,weno_c,weno_zeps,&
                                             evmax,erw,elw,chi,q,fluxes)
   !< Compute convective fluxes along direction `dir`, WENO scheme for space operator.
   integer(I4P), intent(in)    :: dir                                !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                               !< Number of conservative varibales.
   integer(I4P), intent(in)    :: weno_s                             !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                   !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_c(1-weno_s:,1:)               !< Centered polinomials coefficients.
   real(R8P),    intent(in)    :: weno_zeps                          !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)    :: evmax                              !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                      !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                      !< Left  eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: chi                                !< Speed parameter for divergence cleaning.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   integer(I4P)                :: si(3), si_i, si_j, si_k            !< Directional (1=x,2=y,3=z) increment.
   real(R8P)                   :: sir(3)                             !< Directional (1=x,2=y,3=z) increment, real.
   integer(I4P)                :: b, i, j, k                         !< Counter.

   select case(dir)
   case(1)
      si = [1,0,0]
   case(2)
      si = [0,1,0]
   case(3)
      si = [0,0,1]
   endselect
   sir = real(si,R8P)
   si_i = 1-si(1)
   si_j = 1-si(2)
   si_k = 1-si(3)

   do b=1, blocks_number
   do k=si_k, nk
   do j=si_j, nj
   do i=si_i, ni
      call compute_fluxes_convective_ri_weno(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c,                 &
                                             weno_s=weno_s, weno_zeps=weno_zeps,                        &
                                             weno_a=weno_a, weno_p=weno_p, weno_d=weno_d, weno_c=weno_c,&
                                             evmax=evmax,erw=erw,elw=elw,chi=chi,                       &
                                             si=si,sir=sir,q=q,fluxes=fluxes)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_convective_weno

   subroutine compute_fluxes_convective_ri_weno(dir,b,i,j,k,ngc,nv_c,                         &
                                                weno_s,weno_zeps,weno_a,weno_p,weno_d,weno_c, &
                                                evmax,erw,elw,chi,si,sir,q,fluxes)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in)    :: dir                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: b, i, j, k                          !< Counter.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                                !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                              !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_zeps                           !< Parameter to avoid division by zero.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                    !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                 !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                 !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_c(1-weno_s:,1:)                !< Centered polinomials coefficients.
   real(R8P),    intent(in)    :: evmax                               !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                       !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                       !< Left  eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: chi                                 !< Speed parameter for divergence cleaning.
   integer(I4P), intent(in)    :: si(3)                               !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                              !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Fields variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes.
   real(R8P)                   :: fmpc(1:2,1-S_MAX:-1+S_MAX,1:NV_MAX) !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:NV_MAX)                  !< Fluxes +- reconstructed.
   integer(I4P)                :: v, vv                               !< Counter.

   call decompose_fluxes_convective(dir=dir, si=si, sir=sir,                &
                                    b=b, i=i, j=j, k=k, ngc=ngc, nv_c=nv_c, &
                                    weno_s=weno_s, evmax=evmax, elw=elw, chi=chi,    &
                                    q=q, fmpc=fmpc(1:2,1-weno_s:-1+weno_s,1:NV_MAX))
   do v=1, nv_c
      call weno_reconstruct_upwind(S=weno_s, weno_a=weno_a, weno_p=weno_p, weno_d=weno_d,&
                                   weno_zeps=weno_zeps, V=fmpc(1:2,1-weno_s:-1+weno_s,v), VR=fpmr(1:2,v))
   enddo
   ! back projection in conservative variables space
   do v=1, nv_c
      fluxes(v,i,j,k,b) = 0._R8P
      do vv=1,nv_c
         fluxes(v,i,j,k,b) = fluxes(v,i,j,k,b) + erw(vv,v,dir) * (fpmr(1,vv) + fpmr(2,vv))
      enddo
   enddo
   endsubroutine compute_fluxes_convective_ri_weno

   subroutine compute_fluxes_difference(blocks_number, ni, nj, nk, ngc, nv_c, var_Jx, var_Jy, var_Jz, &
                                       dx, dy, dz, flx, fly, flz, q, dq)
   !< Compute fluxes difference.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                            !< Number of conservative varibales in q.
   integer(I4P), intent(in)    :: var_Jx, var_Jy, var_Jz          !< Current variable indices.
   real(R8P),    intent(in)    :: dx(1:), dy(1:), dz(1:)          !< Space steps.
   real(R8P),    intent(in)    :: flx(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)    :: fly(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)    :: flz(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Fields.
   real(R8P),    intent(inout) :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes differences.
   integer(I4P)                :: b, i, j, k, v                   !< Counter.

   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      do v=1, nv_c
         dq(v,i,j,k,b) = - (flx(v,i,j,k,b) - flx(v,i-1,j,k,b)) / dx(b) &
                         - (fly(v,i,j,k,b) - fly(v,i,j-1,k,b)) / dy(b) &
                         - (flz(v,i,j,k,b) - flz(v,i,j,k-1,b)) / dz(b)
      enddo
      ! J sources
      dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) - q(var_Jx,i,j,k,b)
      dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) - q(var_Jy,i,j,k,b)
      dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) - q(var_Jz,i,j,k,b)
      ! corrections
      ! if (d_divergence_cleaner .and. .not.b_divergence_cleaner .and. eta>0._R8P) then
      !    dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
      ! elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
      !    dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
      !    dq(11,i,j,k,b) = dq(11,i,j,k,b) - chi/eta*chi/eta*q(11,i,j,k,b)
      ! endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_difference

   subroutine decompose_fluxes_convective(dir,si,sir,b,i,j,k,ngc,nv_c,weno_s,evmax,elw,q,chi,fmpc)
   !< Decompose convective fluxes.
   integer(I4P), intent(in)    :: dir                           !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: si(3)                         !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                        !< Stencil increment, real cast.
   integer(I4P), intent(in)    :: b, i, j, k                    !< Counter.
   integer(I4P), intent(in)    :: ngc                           !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                          !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                        !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: evmax                         !< Maximum eigenvalue.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                 !< Left eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(in)    :: chi                           !< Speed coefficient for D & B div-cleaning.
   real(R8P),    intent(inout) :: fmpc(1:,1-weno_s:,1:)         !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                   :: fmp(2)                        !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: gc, wc                        !< Increments for fluxes decomposition.
   integer(I4P)                :: v, vv, s, is, js, ks          !< Counter.
   real(R8P)                   :: f(NV_MAX)                     !< Conservative fluxes.
   real(R8P)                   :: p(NV_MAX)                     !< Eigenvalues.

   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      call compute_fluxes_Maxwell(sir=sir,q=q(:,is,js,ks,b),f=f,chi=chi)
      !call compute_eigenvalues_vector(sir=sir, p=p)
      do v=1, nv_c
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv_c
            wc = wc + elw(vv,v,dir) * q(vv,is,js,ks,b)
            gc = gc + elw(vv,v,dir) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax * wc)
         !fmp(2) = 0.5_R8P * (gc + evmax * p(v) * wc)
         fmp(1) = gc - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective

   subroutine impose_div_coil_correction(self, ivar, q)
   !< Impose Constrained Transport Correction on vectorial variable q(ivar:ivar+2).
   !< Note that self%divergence memory is used as buffer, be careful.
   class(prism_cpu_object), intent(inout) :: self                                           !< The equation.
   integer(I4P),            intent(in)    :: ivar                                           !< Variable (start) index in q.
   real(R8P),               intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P)                              :: div_buff(1, 1-self%ngc:self%ni+self%ngc, &
                                                      1-self%ngc:self%nj+self%ngc, &
                                                      1-self%ngc:self%nk+self%ngc, &
                                                      1:self%nb)                            !< Buffer for divergence computation.
   real(R8P)                              :: q_buff(3,1-self%ngc:self%ni+self%ngc, &
                                                      1-self%ngc:self%nj+self%ngc, &
                                                      1-self%ngc:self%nk+self%ngc, &
                                                      1:self%nb)                            !< Buffer for variables computation.
   real(R8P)                              :: grad_buff(3,1-self%ngc:self%ni+self%ngc, &
                                                         1-self%ngc:self%nj+self%ngc, &
                                                         1-self%ngc:self%nk+self%ngc, &
                                                         1:self%nb)                         !< Buffer for gradient computation.
   real(R8P)                              :: dq_buff(1, 1-self%ngc:self%ni+self%ngc, &
                                                     1-self%ngc:self%nj+self%ngc, &
                                                     1-self%ngc:self%nk+self%ngc, &
                                                     1:self%nb)                             !< Buffer for gradient computation.
   real(R8P)                              :: phi_buff(1, 1-self%ngc:self%ni+self%ngc, &
                                                     1-self%ngc:self%nj+self%ngc, &
                                                     1-self%ngc:self%nk+self%ngc, &
                                                     1:self%nb)                             !< Buffer for scalar potential.
   real(R8P)                              :: dq_max                                         !< Maximum residual.
   integer(I4P)                           :: iter                                           !< Counter.
   integer(I4P)                           :: i,j,k,b,v                                      !< Counter.

   !Lascia perdere low memory storage x ora
   !buffer(4,:,:,:,:) è usato come buffer per la correzione di divergenza
   !buffer(5,:,:,:,:) per dq
   !buffer(7,:,:,:,:) per potenziale scalare
   !buffer(4:6,:,:,:,:) è usato come buffer per il gradiente di q

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number)!, buffer=>self%divergence)
   div_buff (:,:,:,:,:) = 0.0_R8P
   grad_buff(:,:,:,:,:) = 0.0_R8P
   dq_buff  (:,:,:,:,:) = 0.0_R8P
   phi_buff (:,:,:,:,:) = 0.0_R8P
   q_buff   (:,:,:,:,:) = q(:,:,:,:,:)

   call self%compute_divergence(ivar=ivar,q=q_buff,divergence=div_buff(1,:,:,:,:))
   print *, ' Max divergence buffered variables: ', maxval(abs(div_buff(:,:,:,:,:)))
   !call self%compute_divergence(ivar=ivar,q=q_buff,divergence=buffer(4,:,:,:,:))
   !print *, ' Max divergence buffered variables 2: ', maxval(abs(buffer(4,:,:,:,:)))
   !buffer(4,:,:,:,:) = div_buff(:,:,:,:)
   !print *, 'Max divergence buffered variables 3: ', maxval(abs(buffer(4,:,:,:,:)))
   if (blocks_number>0) then
      do iter=1, self%flail%iterations
         call compute_smoothing_multigrid(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=1_I4P,blocks_number=blocks_number, &
                                          dxyz=self%field%dxyz,                                           &
                                          f=-div_buff,                                         &
                                          q=phi_buff,                                          &
                                          dq_max=dq_max,                                                  &
                                          dq=dq_buff,                                         &
                                          iterations_init=self%flail%iterations_init,                     &
                                          iterations_fine=self%flail%iterations_fine,                     &
                                          iterations_coarse=self%flail%iterations_coarse)
         call self%compute_gradient(ivar=1,q=phi_buff,gradient=grad_buff)
         !print *, 'valore massima correzione associata all''iterazione ', &
         !         iter, ' - max correction: ', maxval(abs(grad_buff(:,:,:,:,:)))
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do v=1, 3
                        q_buff(v,i,j,k,b) = q_buff(v,i,j,k,b) + grad_buff(v,i,j,k,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         call self%compute_divergence(ivar=ivar,q=q_buff,divergence=div_buff(1,:,:,:,:))
         !print *, 'Coil divergence correction iteration ', iter, ' - max divergence: ', maxval(abs(div_buff(:,:,:,:,:)))
         if (maxval(abs(div_buff(:,:,:,:,:))) < self%flail%tolerance) exit
      enddo
      call self%mpih%print_message('Coil divergence correction converged at iteration '//trim(str(iter,.true.)))
      q(:,:,:,:,:) = q_buff(:,:,:,:,:)
      call self%compute_divergence(ivar=ivar,q=q,divergence=div_buff(1,:,:,:,:))
      print *, 'Coil - max divergence (end Poisson subroutine): ', maxval(abs(div_buff(:,:,:,:,:)))
   endif
   endassociate
   endsubroutine impose_div_coil_correction

	subroutine write_current_behavior_tab(filename, current_density, time)
	character(len=1), parameter  :: TAB = achar(9)
	character(len=*), intent(in) :: filename
	real(R8P),        intent(in) :: current_density
   real(R8P),        intent(in) :: time
	integer(I4P) 					  :: iu, ios

	open(newunit=iu, file=trim(filename), status='unknown', action='write', &
	     form='formatted', position='append', iostat=ios)
	if (ios /= 0) then
	  write(*,'(a,i0)') 'write_current_tab: errore open(), iostat=', ios
	  error stop
	end if
	write(iu,'(ES24.16,a,ES24.16))') time, TAB, current_density
	close(iu)
	endsubroutine write_current_behavior_tab

   subroutine write_single_particle_output(filename, time, q_pic)
	character(len=1), parameter  :: TAB = achar(9)
	character(len=*), intent(in) :: filename
	real(R8P),        intent(in) :: q_pic(1:,1:)
   real(R8P),        intent(in) :: time
	integer(I4P) 					  :: iu, ios, l, j

   l = size(q_pic, dim=1)
  	open(newunit=iu, file=trim(filename), status='unknown', action='write', &
	     form='formatted', position='append', iostat=ios)
	if (ios /= 0) then
	  write(*,'(a,i0)') 'write_current_tab: errore open(), iostat=', ios
	  error stop
	end if
	write(iu,'(ES24.16,8(a,ES24.16))') time, (TAB, q_pic(j,1), j=1,l)
	close(iu)
   endsubroutine

   subroutine compute_field_mean_value(self, q, n_x, n_y, n_z, n_b, mean_value)
   !< Compute mean value of the field in a certain region of the domain.
   class(prism_cpu_object), intent(in)  :: self                                      !< The object
   real(R8P), intent(in)                :: q(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   integer(I4P), intent(in)             :: n_x(2), n_y(2), n_z(2), n_b(2)            !< Number of cells in each direction and number of blocks
   real(R8P), intent(out)               :: mean_value                                !< Mean value of the field out of the fWLayer
   integer(I4P)                         :: i, j, k, b                                !< Counter.

   mean_value = 0._R8P
   do b=n_b(1), n_b(2)
      do k=n_z(1), n_z(2)
         do j=n_y(1), n_y(2)
            do i=n_x(1), n_x(2)
               mean_value = mean_value + q(i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   mean_value = mean_value / real((n_x(2)-n_x(1)+1)*(n_y(2)-n_y(1)+1)*(n_z(2)-n_z(1)+1)*(n_b(2)-n_b(1)+1), R8P)
   endsubroutine compute_field_mean_value

   function dotproduct(a, b) result(dot)
   !< Compute the scalar (dot) product.
   real(R8P), intent(in) :: a(3) !< Left hand side.
   real(R8P), intent(in) :: b(3) !< Left hand side.
   real(R8P)             :: dot  !< Dot product.

   dot = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))
   endfunction dotproduct

   function crossproduct(a, b) result(cross)
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct
endmodule adam_prism_cpu_object

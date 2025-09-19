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

type, extends(prism_common_object) :: prism_cpu_object !commentate procedure AMR e IB
   !< Maxwell equations system class definition, CPU backend.
   real(R8P), allocatable :: flx(:,:,:,:,:) !< Fluxes along x.
   real(R8P), allocatable :: fly(:,:,:,:,:) !< Fluxes along y.
   real(R8P), allocatable :: flz(:,:,:,:,:) !< Fluxes along z.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu !< Allocate CPU data.
      procedure, pass(self) :: initialize   !< Initialize the equation.
      ! IB methods
      procedure, pass(self) :: integrate_eikonal_coils !< Integrate eikonal equation for coils.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_xh5f            !< Save simulation data in XH5F format.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_restart_files   !< Save restart files.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions (and coils) of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      ! numerical methods
      procedure, pass(self) :: compute_dt                 !< Compute time step.
      procedure, pass(self) :: compute_residuals          !< Compute residuals.
      procedure, pass(self) :: compute_residuals_centered !< Compute residuals, centered scheme.
      procedure, pass(self) :: correct_div                !< Correct divergence of q(ivar:2).
      procedure, pass(self) :: integrate                  !< Perform one step integration.
      procedure, pass(self) :: simulate                   !< Perform the simulation.
endtype prism_cpu_object

contains
   ! auxiliary methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.

   call self%mpih%print_message('prism_cpu_object%allocate_cpu start')
   msg_ = self%mpih%myrankstr//'prism_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
              nb=>self%nb, weno_s=>self%weno%S, solids_number=>self%ib%solids_number)
   msg = msg_//' flx '
   call allocate_variable(var=self%flx,ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flx = 0._R8P
   msg = msg_//' fly '
   call allocate_variable(var=self%fly,ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%fly = 0._R8P
   msg = msg_//' flz '
   call allocate_variable(var=self%flz,ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flz = 0._R8P
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
   print '(A)', self%mpih%description()
   call self%mpih%print_message('prism_cpu_object%initialize finish')
   endsubroutine initialize

   ! IB methods
   subroutine integrate_eikonal_coils(self, q)
   !< Integrate eikonal equation.
   class(prism_cpu_object), intent(inout) :: self      !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)     !< Conservative variables.
   integer(I4P)                           :: i_eikonal !< Counter.

   associate(blocks_number=>self%blocks_number, total_coils_number=>self%coil%total_coils_number)
      if (blocks_number > 0) then
         if (total_coils_number > 0) then
            call self%update_ghost(q=q)
            do i_eikonal=1, self%ib%n_eikonal
               call self%mpih%barrier
               call self%ib%evolve_eikonal_coils(q=q, phi=self%coil%phi, n_coils=total_coils_number)
               call self%update_ghost(q=q)
            enddo
            !call self%ib%invert_eikonal_coils(q=q)
            call self%mpih%barrier
         endif
      endif
   endassociate
   endsubroutine integrate_eikonal_coils

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save simulation data in HDF5 format.
   class(prism_cpu_object), intent(inout)        :: self             !< The equation.
   character(*),            intent(in), optional :: output_basename  !< Output basename.
   logical,                 intent(in), optional :: with_ghost       !< Flag to save ghost cells.
   character(:), allocatable                     :: output_basename_ !< Output basename, local var.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it,9))
   if (present(output_basename)) output_basename_ = trim(output_basename)
   call self%adam%io%save_xh5f(basename=trim(output_basename_), &
                               q=self%q, q_name=self%q_name,    &
                               with_ghost=with_ghost,           &
                               with_cell_morton=.true.,         &
                               t=self%time%it, time=self%time%time)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_xh5f

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

   subroutine save_restart_files(self)
   !< Save restart files.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q=self%q)

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max)) then
         if (.not.self%physics%D_divergence_cleaner .and. .not.self%physics%B_Divergence_cleaner) then
            call self%slices%save_mat(basename=self%io%output_basename, &
                                      it=self%time%it,                  &
                                      it_max=self%time%it_max,          &
                                      time=self%time%time,              &
                                      time_max=self%time%time_max,      &
                                      adam=self%adam,                   &
                                      q=self%q,                         &
                                      q_name=['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz '])
         elseif (self%physics%div_corr_var == DIV_CORR_VAR_HYPER .and. self%physics%D_divergence_cleaner & 
                  .and. .not.self%physics%B_Divergence_cleaner) then
            call self%slices%save_mat(basename=self%io%output_basename, &
                                      it=self%time%it,                  &
                                      it_max=self%time%it_max,          &
                                      time=self%time%time,              &
                                      time_max=self%time%time_max,      &
                                      adam=self%adam,                   &
                                      q=self%q,                         &
                                      q_name=['Dx  ','Dy  ','Dz  ','Bx  ','By  ','Bz  ','Jx  ','Jy  ','Jz  ','Phi '])
         elseif (self%physics%div_corr_var == DIV_CORR_VAR_HYPER .and. self%physics%D_divergence_cleaner &
                  .and. self%physics%B_Divergence_cleaner) then
            call self%slices%save_mat(basename=self%io%output_basename, &
                                      it=self%time%it,                  &
                                      it_max=self%time%it_max,          &
                                      time=self%time%time,              &
                                      time_max=self%time%time_max,      &
                                      adam=self%adam,                   &
                                      q=self%q,                         &
                                      q_name=['Dx  ','Dy  ','Dz  ','Bx  ','By  ','Bz  ','Jx  ','Jy  ','Jz  ','Phi ','Psi '])
         endif
      endif
   endif
   endsubroutine save_simulation_data

   ! IC/BC
   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(prism_cpu_object), intent(in)    :: self                 !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:)    !< Conservative variables.
   integer(I4P)                        :: b, c, i, j, k, v        !< Counter.
   integer(I4P)                        :: idelta,jdelta,kdelta    !< IJK delta step for extrapolation.
   integer(I4P)                        :: bc_type                 !< Boundary condition type.
   integer(I4P)                        :: crown                   !< Crown counter.
   integer(I4P)                        :: fec                     !< Boundary fec (1 to 26).
   integer(I4P)                        :: fec_1_6                 !< Boundary fec (1 to 6).
   integer(I4P)                        :: alfa_D, beta_D, gamma_D !< Indici alfa beta gamma come in Barbas.
   integer(I4P)                        :: alfa_B, beta_B, gamma_B !< Indici alfa beta gamma come in Barbas.
   real(R8P)                           :: s1                      !< Coefficiente pari a +-1.
   real(R8P)                           :: ds                      !< Distanza tra le celle in x, y o z.
   real(R8P)                           :: ngc_r, crown_r          !< Numero di gc totale, reale
   real(R8P)                           :: ref(1:9)                !< Vettore di stato di riferimento per assegnazione gc.
   real(R8P)                           :: fi, f                   !< Variabili phi e f fWL.

   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown, &
             nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), &
             dz=>self%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             D_divergence_cleaner=>self%physics%D_divergence_cleaner, dt=>self%time%dt, &
             B_divergence_cleaner=>self%physics%B_divergence_cleaner, chi=>self%physics%chi)
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
                  do v=1, 9
                     q(v,i,j,k,b) = q(v,i-idelta,j-jdelta,k-kdelta,b) !ni,j,k coordinate della cella da cui prendo i valori
                  enddo
                  if (self%physics%D_divergence_cleaner) then
                     q(10,i,j,k,b) = 0._R8P
                     !q(10,i,j,k,b) = q(10,i-idelta,j-jdelta,k-kdelta,b) - dx(b)/((chi*sqrt(1/(MU0*EPS0)))*dt)* &
                                    !(q(10,i-idelta,j-jdelta,k-kdelta,b)-q_old(10,i-idelta,j-jdelta,k-kdelta,b))
                  endif
                  if (self%physics%B_divergence_cleaner) then
                     q(11,i,j,k,b) = 0._R8P
                     !q(11,i,j,k,b) = q(11,i-idelta,j-jdelta,k-kdelta,b) - dx(b)/((chi*sqrt(1/(MU0*EPS0)))*dt)* &
                                    !(q(11,i-idelta,j-jdelta,k-kdelta,b)-q_old(11,i-idelta,j-jdelta,k-kdelta,b))
                  endif
               elseif (bc_type == BC_fWLayer) then
                  !print *, fec
                  if (fec <= 6) then
                     select case(fec)
                     !Identifico gli alfa beta gamma come nel paper di Barbas, distinguendo tra alfa_D e alfa_B ecc

                     case(1) ! x- face alfa = 2, beta = 3, gamma = 1
                        s1 = 1.0_R8P
                        alfa_D = 2_I4P
                        beta_D = 3_I4P
                        gamma_D = 1_I4P
                        alfa_B = 5_I4P
                        beta_B = 6_I4P
                        gamma_B = 4_I4P
                        ds = dx(b) !distanza tra le celle in x
                        ref = q(:,1,j,k,b) !vettore di stato di riferimento per assegnazione gc

                     case(2) ! x+ face
                        s1 = -1.0_R8P
                        alfa_D = 2_I4P
                        beta_D = 3_I4P
                        gamma_D = 1_I4P
                        alfa_B = 5_I4P
                        beta_B = 6_I4P
                        gamma_B = 4_I4P
                        ds = dx(b) !distanza tra le celle in x
                        ref = q(:,ni,j,k,b)

                     case(3) ! y- face
                        s1 = 1.0_R8P
                        alfa_D = 3_I4P
                        beta_D = 1_I4P
                        gamma_D = 2_I4P
                        alfa_B = 6_I4P
                        beta_B = 4_I4P
                        gamma_B = 5_I4P
                        ds = dy(b) !distanza tra le celle in y
                        ref = q(:,i,1,k,b)

                     case(4) ! y+ face
                        s1 = -1.0_R8P
                        alfa_D = 3_I4P
                        beta_D = 1_I4P
                        gamma_D = 2_I4P
                        alfa_B = 6_I4P
                        beta_B = 4_I4P
                        gamma_B = 5_I4P
                        ds = dy(b) !distanza tra le celle in y
                        ref = q(:,i,nj,k,b)

                     case(5) ! z- face
                        s1 = 1.0_R8P
                        alfa_D = 1_I4P
                        beta_D = 2_I4P
                        gamma_D = 3_I4P
                        alfa_B = 4_I4P
                        beta_B = 5_I4P
                        gamma_B = 6_I4P
                        ds = dz(b) !distanza tra le celle in z
                        ref = q(:,i,j,1,b)

                     case(6) ! z+ face
                        s1 = -1.0_R8P
                        alfa_D = 1_I4P
                        beta_D = 2_I4P
                        gamma_D = 3_I4P
                        alfa_B = 4_I4P
                        beta_B = 5_I4P
                        gamma_B = 6_I4P
                        ds = dz(b) !distanza tra le celle in z
                        ref = q(:,i,j,nk,b)

                     endselect
                     ngc_r = real(ngc,R8P)
                     crown_r = real(crown, R8P)
                     if (ngc < 40_I4P) then
                        fi = 1/150._R8P*(-7.0_R8P*ngc_r**2 + 255._R8P*ngc_r + 250._R8P) !polinomio di Barbas
                     else
                        fi = 25.0_R8P
                     endif
                     !x - xa è la distanza tra il centro della gc considerata e il bordo del dominio (fatto col centro cella, vedremo)
                     !è pari quindi a (ngc_r - crown_r) * ds

                     !xb - xa è la distanza tra il centro della gc più esterna considerata e il bordo del dominio (fatto col centro cella, vedremo)
                     !è pari quindi a C * ds

                     f = 1._R8P/fi*LOG10(((ngc_r-crown_r)*ds)/(ngc_r*ds)*(10._R8P**fi-1._R8P)+1._R8P) !funzione f
      
                     q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s1*(f-1._R8P)*ref(beta_B)*EPS0**0.5_R8P + &
                                          (f+1._R8P)*ref(alfa_D)*MU0**0.5_R8P)

                     q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f+1._R8P)*ref(beta_B)*EPS0**0.5_R8P + &
                                          s1*(f-1._R8P)*ref(alfa_D)*MU0**0.5_R8P)

                     q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s1*(f-1._R8P)*ref(alfa_B)*EPS0**0.5_R8P + &
                                          (f+1._R8P)*ref(beta_D)*MU0**0.5_R8P)

                     q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f+1._R8P)*ref(alfa_B)*EPS0**0.5_R8P - &
                                          s1*(f-1._R8P)*ref(beta_D)*MU0**0.5_R8P)

                     q(gamma_D,i,j,k,b) = ref(gamma_D)

                     q(gamma_B,i,j,k,b) = ref(gamma_B)

                     q(7:9,i,j,k,b) = 0._R8P
                  else
                     do v=1, nv
                       q(v,i,j,k,b) = 0.0_R8P
                     enddo
                  endif
               elseif (bc_type == BC_Silver_Muller) then
                  !print *, fec
                  if (fec <= 6) then
                     select case(fec)
                     !Identifico gli alfa beta gamma come nel paper di Barbas, distinguendo tra alfa_D e alfa_B ecc

                     case(1) ! x- face alfa = 2, beta = 3, gamma = 1
                        s1 = 1.0_R8P
                        alfa_D = 2_I4P
                        beta_D = 3_I4P
                        gamma_D = 1_I4P
                        alfa_B = 5_I4P
                        beta_B = 6_I4P
                        gamma_B = 4_I4P
                        ds = dx(b) !distanza tra le celle in x
                        ref = q(:,i+1,j,k,b) !vettore di stato di riferimento per assegnazione gc

                     case(2) ! x+ face
                        s1 = -1.0_R8P
                        alfa_D = 2_I4P
                        beta_D = 3_I4P
                        gamma_D = 1_I4P
                        alfa_B = 5_I4P
                        beta_B = 6_I4P
                        gamma_B = 4_I4P
                        ds = dx(b) !distanza tra le celle in x
                        ref = q(:,i-1,j,k,b)

                     case(3) ! y- face
                        s1 = 1.0_R8P
                        alfa_D = 3_I4P
                        beta_D = 1_I4P
                        gamma_D = 2_I4P
                        alfa_B = 6_I4P
                        beta_B = 4_I4P
                        gamma_B = 5_I4P
                        ds = dy(b) !distanza tra le celle in y
                        ref = q(:,i,j+1,k,b)

                     case(4) ! y+ face
                        s1 = -1.0_R8P
                        alfa_D = 3_I4P
                        beta_D = 1_I4P
                        gamma_D = 2_I4P
                        alfa_B = 6_I4P
                        beta_B = 4_I4P
                        gamma_B = 5_I4P
                        ds = dy(b) !distanza tra le celle in y
                        ref = q(:,i,j-1,k,b)

                     case(5) ! z- face
                        s1 = 1.0_R8P
                        alfa_D = 1_I4P
                        beta_D = 2_I4P
                        gamma_D = 3_I4P
                        alfa_B = 4_I4P
                        beta_B = 5_I4P
                        gamma_B = 6_I4P
                        ds = dz(b) !distanza tra le celle in z
                        ref = q(:,i,j,k+1,b)

                     case(6) ! z+ face
                        s1 = -1.0_R8P
                        alfa_D = 1_I4P
                        beta_D = 2_I4P
                        gamma_D = 3_I4P
                        alfa_B = 4_I4P
                        beta_B = 5_I4P
                        gamma_B = 6_I4P
                        ds = dz(b) !distanza tra le celle in z
                        ref = q(:,i,j,k-1,b)

                     endselect
                     ngc_r = real(ngc,R8P)
                     crown_r = real(crown, R8P)

                     ! fWLayer con f = 0 è Silver-Muller
                     f = 0._R8P !funzione f

                     q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s1*(f-1._R8P)*ref(beta_B)*EPS0**0.5_R8P + &
                                          (f+1._R8P)*ref(alfa_D)*MU0**0.5_R8P)

                     q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f+1._R8P)*ref(beta_B)*EPS0**0.5_R8P + &
                                          s1*(f-1._R8P)*ref(alfa_D)*MU0**0.5_R8P)

                     q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s1*(f-1._R8P)*ref(alfa_B)*EPS0**0.5_R8P + &
                                          (f+1._R8P)*ref(beta_D)*MU0**0.5_R8P)

                     q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f+1._R8P)*ref(alfa_B)*EPS0**0.5_R8P - &
                                          s1*(f-1._R8P)*ref(beta_D)*MU0**0.5_R8P)

                     q(gamma_D,i,j,k,b) = ref(gamma_D)

                     q(gamma_B,i,j,k,b) = ref(gamma_B)

                     q(7:9,i,j,k,b) = 0._R8P
                  else
                     do v=1, nv
                        q(v,i,j,k,b) = 0.0_R8P
                     enddo
                  endif
               elseif (bc_type == BC_EXTRAP_DIRICHLET) then
                     do v=1, 9
                        q(v,i,j,k,b) = 0._R8P
                     enddo
               endif

            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions and coils on field.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field, q=self%q)
   call self%coil%set_coils(physics=self%physics, field=self%field)
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(prism_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),               intent(inout)        :: q(1:,         &
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1:)           !< Conservative variables.
   integer(I4P),            intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
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
   if (do_set_bc)       call self%set_boundary_conditions(q=q)
   endsubroutine update_ghost

   ! numerical methods
   subroutine compute_dt(self)
   class(prism_cpu_object), intent(inout) :: self                            !< The equation.
   real(R8P)                              :: umax                            !< Maximum speed of waves propagation (light speed).
   real(R8P)                              :: dxyz_min                        !< Minimal space step.
   real(R8P)                              :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                           :: b                               !< Counter.

   dxyz_min = huge(1._R8P)
   associate(blocks_number=>self%blocks_number, dxyz=>self%field%dxyz, d_divergence_cleaner=>self%physics%d_divergence_cleaner,&
             chi=>self%physics%chi, eta=>self%physics%eta)
   call compute_dxyz_min(blocks_number=blocks_number, dxyz=dxyz, dxyz_min=dxyz_min)
   umax = C0
   if (d_divergence_cleaner .and. self%physics%div_corr_var == DIV_CORR_VAR_HYPER) umax = max(chi*C0, eta*C0)
   self%time%dt = self%time%CFL*dxyz_min / umax
   endassociate
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endsubroutine compute_dt

   subroutine compute_residuals(self, q, dq)
   !< Compute residuals of equation.
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

   call self%update_ghost(q=q)
   !call self%integrate_eikonal_coils(q=q)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:),                         &
             flx=>self%flx, fly=>self%fly, flz=>self%flz,                                                          &
             weno_s=>self%weno%S, weno_zeps=>self%weno%zeps,                                                       &
             weno_a=>self%weno%a, weno_p=>self%weno%p, weno_d=>self%weno%d,                                        &
             evmax=>self%physics%evmax, erw=>self%physics%erw, elw=>self%physics%elw)
   if (blocks_number > 0) then
      call compute_fluxes_convective(dir=1,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,      &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,                                                &
                                     q=q,fluxes=flx)
      call compute_fluxes_convective(dir=2,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,      &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,                                                &
                                     q=q,fluxes=fly)
      call compute_fluxes_convective(dir=3,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,      &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,                                                &
                                     q=q,fluxes=flz)
      call compute_fluxes_difference(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv_c=nv_c, &
                                     dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, dq=dq, q=q)
   endif
   endassociate
   endsubroutine compute_residuals

   subroutine compute_residuals_centered(self, q, dq)
   !< Compute residuals of equation, centerd scheme.
   class(prism_cpu_object), intent(inout) :: self                 !< The equation.
   real(R8P),               intent(inout) :: q(1:,       &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)                !< Conservative variables.
   real(R8P),               intent(inout) :: dq(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:)               !< Residuals.
   integer(I4P)                           :: b, i, j, k, v        !< Counter.

   call self%update_ghost(q=q)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:),                         &
             flx=>self%flx, fly=>self%fly, flz=>self%flz)
   
   

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               call compute_convective_fluxes_maxwell(sir=[1.0_R8P, 0.0_R8P, 0.0_R8P], q=q(:,i,j,k,b), f=flx(:,i,j,k,b))
               call compute_convective_fluxes_maxwell(sir=[0.0_R8P, 1.0_R8P, 0.0_R8P], q=q(:,i,j,k,b), f=fly(:,i,j,k,b))
               call compute_convective_fluxes_maxwell(sir=[0.0_R8P, 0.0_R8P, 1.0_R8P], q=q(:,i,j,k,b), f=flz(:,i,j,k,b))
            enddo
         enddo
      enddo
   enddo

   do b=1, blocks_number
      do v=1, nv_c
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  dq(v,i,j,k,b) = -( (flx(v,i+1,j,k,b)-flx(v,i-1,j,k,b))/(2*dx(b)) + &
                                     (fly(v,i,j+1,k,b)-fly(v,i,j-1,k,b))/(2*dy(b)) + &
                                     (flz(v,i,j,k+1,b)-flz(v,i,j,k-1,b))/(2*dz(b)))
               enddo
            enddo
         enddo
      enddo
   enddo

   !J sources
   dq(VAR_DX,:,:,:,:) = dq(VAR_DX,:,:,:,:) - q(VAR_JX,:,:,:,:)
   dq(VAR_DY,:,:,:,:) = dq(VAR_DY,:,:,:,:) - q(VAR_JY,:,:,:,:)
   dq(VAR_DZ,:,:,:,:) = dq(VAR_DZ,:,:,:,:) - q(VAR_JZ,:,:,:,:)

   endassociate
   endsubroutine compute_residuals_centered

   subroutine integrate(self, do_ghost_syncro)
   !< Perform one step integration.
   class(prism_cpu_object), intent(inout)         :: self             !< The equation.
   logical,                 intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   logical                                        :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number,   &
             time=>self%time%time, A=>self%coil%A, f=>self%coil%f, phase=>self%coil%phase,              &
             coil_flag =>self%coil%coil_flag, d=>self%coil%d, td=>self%coil%td, j_vec=>self%coil%j_vec, &
             dx=>self%field%dxyz(1,1), dxyz=>self%field%dxyz)

   if (self%coil%total_coils_number >= 1_I4P) then
      call compute_coils_current(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, time=time, A=A, d=d, &
                                 f=f, phase=phase, coil_flag=coil_flag, td=td, j_vec=j_vec, dx=dx, q=self%q)
   endif

   call self%rk%initialize_stages(q=self%q)

   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      ! low storage RK working on q_rk(:,:,:,:,:,1)/q_gpu as stages, update q_gpu in place
      do s=1, self%rk%nrk

         call self%compute_residuals_centered(q=self%q, dq=self%dq)

         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk%compute_stage_ls(s=s,dt=self%time%dt,phi=self%ib%phi,dq=self%dq,q=self%q)
         else
            call self%rk%compute_stage_ls(s=s,dt=self%time%dt,dq=self%dq,q=self%q)
         endif

      enddo
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      ! RK working on q_rk as stages
      do s=1, self%rk%nrk
         if (self%ib%solids_number>0) then
            call self%rk%compute_stage(s=s, dt=self%time%dt, phi=self%ib%phi)
         else
            call self%rk%compute_stage(s=s, dt=self%time%dt)
         endif
         call self%compute_residuals_centered(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk%assign_stage(s=s, q=self%dq, phi=self%ib%phi)
         else
            call self%rk%assign_stage(s=s, q=self%dq)
         endif
      enddo
      if (self%ib%solids_number>0) then
         call self%rk%update_q(dt=self%time%dt, phi=self%ib%phi, q=self%q)
      else
         call self%rk%update_q(dt=self%time%dt, q=self%q)
      endif
   endselect
   if (self%physics%div_corr_var == 'POISSON' .and. self%physics%D_divergence_cleaner) then
      call self%correct_div(ivar=1_I4P) ! correct div(D)
   endif
      if (self%physics%div_corr_var == 'POISSON' .and. self%physics%D_divergence_cleaner) then
      call self%correct_div(ivar=1_I4P) ! correct div(D)
   endif
   call compute_div(ni=ni,nj=nj,nk=nk,ngc=ngc,blocks_number=blocks_number,dxyz=dxyz,ivar=1,q=self%q,div=self%field_div(1,:,:,:,:))
   call compute_div(ni=ni,nj=nj,nk=nk,ngc=ngc,blocks_number=blocks_number,dxyz=dxyz,ivar=4,q=self%q,div=self%field_div(2,:,:,:,:))
   call compute_div(ni=ni,nj=nj,nk=nk,ngc=ngc,blocks_number=blocks_number,dxyz=dxyz,ivar=7,q=self%q,div=self%field_div(3,:,:,:,:))
   endassociate
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(prism_cpu_object), intent(inout) :: self             !< The equation.
   character(*),            intent(in)    :: filename         !< Input file name.
   real(R8P)                              :: timing(1:2)      !< Tic toc timing.
   real(R8P)                              :: timing_step(1:2) !< Tic toc timing.
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
   call self%save_simulation_data

   if (self%mpih%myrank==0) call self%io%open_file_residuals(nv=self%nv)

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

      if (((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
         ((self%time%it>=self%time%it_max).and.(self%time%it_max > 0))) exit integration

      call self%mpih%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
   enddo integration
   call self%mpih%barrier(tictoc=.true., timing=timing(2), single=.true.)
   call self%save_simulation_data
   if (self%mpih%myrank==0) call self%io%close_file_residuals
   call self%mpih%finalize
   endsubroutine simulate

   subroutine correct_div(self, ivar)
   !< Correct divergence of q(ivar:2).
   class(prism_cpu_object), intent(inout) :: self                       !< The equation.
   integer(I4P),            intent(in)    :: ivar                       !< Variable (start) index in q.
   real(R8P)                              :: div(1:1,                        &
                                                 1-self%ngc:self%ni+self%ngc,&
                                                 1-self%ngc:self%nj+self%ngc,&
                                                 1-self%ngc:self%nk+self%ngc,&
                                                 1:self%blocks_number)  !< Divergence of q(ivar).
   real(R8P)                              ::  dq(1:1,                        &
                                                 1-self%ngc:self%ni+self%ngc,&
                                                 1-self%ngc:self%nj+self%ngc,&
                                                 1-self%ngc:self%nk+self%ngc,&
                                                 1:self%blocks_number)  !< Residual of smoothing.
   real(R8P)                              :: grad(1:3,                        &
                                                  1-self%ngc:self%ni+self%ngc,&
                                                  1-self%ngc:self%nj+self%ngc,&
                                                  1-self%ngc:self%nk+self%ngc,&
                                                  1:self%blocks_number) !< Gradient of phi.
   real(R8P)                              :: dq_max                     !< Maximum residual.
   integer(I4P)                           :: iter                       !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number, dxyz=>self%field%dxyz)
   call compute_div(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, dxyz=dxyz, ivar=ivar, &
                    q=self%q, div=div(1,:,:,:,:))
   if (blocks_number>0) then
      do iter=1, self%flail%iterations
         call compute_smoothing_multigrid(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=1_I4P,blocks_number=blocks_number, &
                                          dxyz=dxyz,                                                      &
                                          f=-div(1:1,:,:,:,:),                                            &
                                          q=self%phid,                                                    &
                                          dq_max=dq_max,                                                  &
                                          dq=dq(1:1,:,:,:,:),                                             &
                                          iterations_init=self%flail%iterations_init,                     &
                                          iterations_fine=self%flail%iterations_fine,                     &
                                          iterations_coarse=self%flail%iterations_coarse)
         if (dq_max < self%flail%tolerance) exit
      enddo
      call self%mpih%print_message('FLAIL convergence reached at iteration '//trim(str(iter,.true.)))
      call compute_grad(ni=ni,nj=nj,nk=nk,ngc=ngc,blocks_number=blocks_number,dxyz=dxyz,ivar=1,q=self%phid,grad=grad)
      self%q(ivar:ivar+2,1:ni,1:nj,1:nk,1:blocks_number) = self%q(ivar:ivar+2,1:ni,1:nj,1:nk,1:blocks_number) &
                                                         + grad(  1:3,        1:ni,1:nj,1:nk,1:blocks_number)
   endif
   endassociate
   endsubroutine correct_div

   ! non TBP
   subroutine compute_div(ni, nj, nk, ngc, blocks_number, dxyz, ivar, q, div)
   !< Compute div(q(ivar). Finite difference central scheme.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                     !< Space steps.
   integer(I4P), intent(in)    :: ivar                            !< Variable (vectorial) of q.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Field variables.
   real(R8P),    intent(inout) :: div(   1-ngc:,1-ngc:,1-ngc:,1:) !< Divergence of D, B.
   integer(I4P)                :: i,j,k,b                         !< Counter

   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      div(i,j,k,b) = 0.5_R8P*((q(ivar  ,i+1,j,k,b) - q(ivar  ,i-1,j,k,b))/dxyz(1,b) + &
                              (q(ivar+1,i,j+1,k,b) - q(ivar+1,i,j-1,k,b))/dxyz(2,b) + &
                              (q(ivar+2,i,j,k+1,b) - q(ivar+2,i,j,k-1,b))/dxyz(3,b))

   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_div

   subroutine compute_grad(ni, nj, nk, ngc, blocks_number, dxyz, q, ivar, grad)
   !< Compute gradient of q(ivar). Finite difference central scheme.
   integer(I4P), intent(in)    :: ni                               !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                               !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                               !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                              !< Ghost cells number.
   integer(I4P), intent(in)    :: blocks_number                    !< Number of current blocks.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                      !< Space steps.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Field component to which apply gradient.
   integer(I4P), intent(in)    :: ivar                             !< Index of variable for computing the gradient.
   real(R8P),    intent(inout) :: grad(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Gradient of q(ivar).
   integer(I4P)                :: i, j, k, b                       !< Counter.

   !$omp parallel do collapse(4) default(firstprivate) shared(q,grad)
   do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         grad(1,i,j,k,b) = (q(ivar,i+1,j,k,b) - q(ivar,i-1,j,k,b))/(2*dxyz(1,b))
         grad(2,i,j,k,b) = (q(ivar,i,j+1,k,b) - q(ivar,i,j-1,k,b))/(2*dxyz(2,b))
         grad(3,i,j,k,b) = (q(ivar,i,j,k+1,b) - q(ivar,i,j,k-1,b))/(2*dxyz(3,b))
      enddo
      enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_grad

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

   subroutine compute_fluxes_convective(dir,blocks_number,ni,nj,nk,ngc,nv_c,weno_s,weno_a,weno_p,weno_d,weno_zeps,&
                                        evmax,erw,elw,q,fluxes)
   !< Compute convective fluxes along direction `dir`.
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
   real(R8P),    intent(in)    :: weno_zeps                          !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)    :: evmax                              !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                      !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                      !< Left  eigenvectors for WENO reconstruction.
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
      call compute_fluxes_convective_ri(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c,  &
                                        weno_s=weno_s, weno_zeps=weno_zeps,         &
                                        weno_a=weno_a, weno_p=weno_p, weno_d=weno_d,&
                                        evmax=evmax,erw=erw,elw=elw,                &
                                        si=si,sir=sir,q=q,fluxes=fluxes)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_convective

   subroutine compute_fluxes_convective_ri(dir,b,i,j,k,ngc,nv_c,                  &
                                           weno_s,weno_zeps,weno_a,weno_p,weno_d, &
                                           evmax,erw,elw,si,sir,q,fluxes)
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
   real(R8P),    intent(in)    :: evmax                              !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                      !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                      !< Left  eigenvectors for WENO reconstruction.
   integer(I4P), intent(in)    :: si(3)                               !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                              !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Fields variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes.
   real(R8P)                   :: fmpc(1:2,1-S_MAX:-1+S_MAX,1:NV_MAX) !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:NV_MAX)                  !< Fluxes +- reconstructed.
   integer(I4P)                :: v, vv                               !< Counter.

   call decompose_fluxes_convective(dir=dir, si=si, sir=sir,                &
                                    b=b, i=i, j=j, k=k, ngc=ngc, nv_c=nv_c, &
                                    weno_s=weno_s, evmax=evmax, elw=elw,    &
                                    q=q, fmpc=fmpc)
   do v=1, nv_c
      call weno_reconstruct_upwind(S=weno_s, weno_a=weno_a, weno_p=weno_p, weno_d=weno_d,&
                                   weno_zeps=weno_zeps, V=fmpc(:,:,v), VR=fpmr(:,v))
   enddo
   ! back projection in conservative variables space
   do v=1, nv_c
      fluxes(v,i,j,k,b) = 0._R8P
      do vv=1,nv_c
         fluxes(v,i,j,k,b) = fluxes(v,i,j,k,b) + erw(vv,v,dir) * (fpmr(1,vv) + fpmr(2,vv))
      enddo
   enddo
   endsubroutine compute_fluxes_convective_ri

   subroutine compute_fluxes_difference(blocks_number, ni, nj, nk, ngc, nv_c, dx, dy, dz, flx, fly, flz, q, dq)
   !< Compute fluxes difference.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                            !< Number of conservative varibales in q.
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
      dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) - q(VAR_JX,i,j,k,b)
      dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) - q(VAR_JY,i,j,k,b)
      dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) - q(VAR_JZ,i,j,k,b)
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

   subroutine decompose_fluxes_convective(dir,si,sir,b,i,j,k,ngc,nv_c,weno_s,evmax,elw,q,fmpc)
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
   real(R8P),    intent(inout) :: fmpc(1:,1-weno_s:,1:)         !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                   :: fmp(2)                        !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: gc, wc                        !< Increments for fluxes decomposition.
   integer(I4P)                :: v, vv, s, is, js, ks          !< Counter.
   real(R8P)                   :: f(NV_MAX)                     !< Conservative fluxes.

   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      call compute_convective_fluxes_maxwell(sir=sir,q=q(:,is,js,ks,b),f=f)
      do v=1, nv_c
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv_c
            wc = wc + elw(vv,v,dir) * q(vv,is,js,ks,b)
            gc = gc + elw(vv,v,dir) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax * wc)
         fmp(1) = gc - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective

   subroutine compute_coils_current(ni, nj, nk, ngc, blocks_number, time, A, d, f, phase, coil_flag, td, j_vec, dx, q)
   !< Compute current coils sources.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   real(R8P),    intent(in)    :: time                               !< Simulation time, to compute current value if AC.
   real(R8P),    intent(in)    :: A(0:)                              !< Current amplitude (A).
   real(R8P),    intent(in)    :: d(0:)                              !< Wire diameter.
   real(R8P),    intent(in)    :: f(0:)                              !< Current frequency, if AC (Hz).
   real(R8P),    intent(in)    :: phase(0:)                          !< Current initial phase, if AC.
   integer(I4P), intent(in)    :: coil_flag(1-ngc:,1-ngc:,1-ngc:,1:) !< Coils ID map.
   real(R8P),    intent(in)    :: td                                 !< Coils transitory delay.
   real(R8P),    intent(in)    :: j_vec(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Current J versors into coils.
   real(R8P),    intent(in)    :: dx                                 !< Space step in x direction.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   real(R8P)                   :: current_density                    !< Current density.
   real(R8P)                   :: g                                  !< Starting polynomial transitory of coils.
   integer(I4P)                :: w_, w_c_                           !< Step function coeff to avoid if in parallel regions.
   real(R8P)                   :: g_, f_                             !< Current coefficients.
   integer(I4P)                :: coil_id                            !< Uniq coild ID.
   integer(I4P)                :: i,j,k,b                            !< Counter.

   !if (time >= td) then
   !   q(VAR_JX,:,:,:,:) = 0._R8P
   !   q(VAR_JY,:,:,:,:) = 0._R8P
   !   q(VAR_JZ,:,:,:,:) = 0._R8P
   !else  
      g = 10._R8P*(time/td)**3 - 15._R8P*(time/td)**4 + 6._R8P*(time/td)**5
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         coil_id = coil_flag(i,j,k,b)

         ! use step function to avoid the following original if
         !if (time < td) then
         !   current_density = g*A(coil_id)/((d(coil_id)-dx)**2)*cos(phase(coil_id)*pi/180.0_R8P)
         !else
         !   current_density = A(coil_id)/((d(coil_id)-dx)**2)*cos(2*pi*f(coil_id)*(time-td) + &
         !   phase(coil_id)*pi/180.0_R8P)
         !endif
         w_   = nint(sign(1._R8P,td-time) + 1._R8P)/2   ! = 1 if td>time, = 0                            if td<time
         w_c_ = 1_I4P - w_                            ! = 0 if td>time, = 1                              if td<time
         g_   = w_ * g + w_c_                         ! = g if td>time, = 1                              if td<time
         f_   = w_c_ * 2._R8P*PI*f(coil_id)*(time-td) ! = 0 if td>time, = 2._R8P*PI*f(coil_id)*(time-td) if td<time
         !current_density = g_ * A(coil_id) / ((d(coil_id))**2) * cos(f_ + phase(coil_id)*PI/180.0_R8P)
         current_density = g_ * A(coil_id) * cos(f_ + phase(coil_id)*PI/180.0_R8P)*j_vec(4,i,j,k,b)
         !if (coil_id == 1_I4P) then
         !   print*, A(coil_id)
         !   print*, current_density
         !   print*, j_vec(4,i,j,k,b)
         !endif
         ! the following if is not necessary because j_vec is zero everywhere except in coils
         if (coil_id /= 0_I4P) then
            q(VAR_JX,i,j,k,b) = current_density * j_vec(1,i,j,k,b)
            q(VAR_JY,i,j,k,b) = current_density * j_vec(2,i,j,k,b)
            q(VAR_JZ,i,j,k,b) = current_density * j_vec(3,i,j,k,b)
            !print*, q(VAR_JX,i,j,k,b), q(VAR_JY,i,j,k,b), q(VAR_JZ,i,j,k,b)
         endif
      enddo
      enddo
      enddo
      enddo
   !endif
   endsubroutine compute_coils_current

   function sq_norm(a) result(sq)
   !< Return the square of the norm of vector.
   real(R8P), intent(in)  :: a(3)     !< Input vector
   real(R8P)              :: sq       !< Square norm of input

   sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
   endfunction sq_norm
endmodule adam_prism_cpu_object

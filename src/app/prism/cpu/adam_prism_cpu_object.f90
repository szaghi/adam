!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, CPU backend.
module adam_prism_cpu_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, CPU backend.

use adam_common_library
use adam_prism_common_library
use penf
use mpi

implicit none
private
public :: prism_cpu_object

type, extends(prism_common_object) :: prism_cpu_object !commentate procedure AMR e IB
   !< Maxwell equations system class definition, CPU backend.
   real(R8P), allocatable :: dq( :,:,:,:,:) !< Residuals right hand side.
   real(R8P), allocatable :: flx(:,:,:,:,:) !< Fluxes along x.
   real(R8P), allocatable :: fly(:,:,:,:,:) !< Fluxes along y.
   real(R8P), allocatable :: flz(:,:,:,:,:) !< Fluxes along z.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu !< Allocate CPU data.
      procedure, pass(self) :: initialize   !< Initialize the equation.
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
      procedure, pass(self) :: compute_dt        !< Compute time step.
      procedure, pass(self) :: compute_residuals !< Compute residuals.
      procedure, pass(self) :: integrate         !< Perform one step integration.
      procedure, pass(self) :: simulate          !< Perform the simulation.
endtype prism_cpu_object

interface assign_omp
!< Assign array to scalar value with OpenMP threads.
module procedure assign_omp_R8P_5D
endinterface assign_omp

contains
   ! auxiliary methods
   subroutine allocate_cpu(self) !modificata associazione eliminando ns e n_aux
   !< Allocate CPU data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.

   call self%mpih%print_message('prism_cpu_object%allocate_cpu start')
   msg_ = self%mpih%myrankstr//'prism_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
              nb=>self%nb, weno_s=>self%weno%S, solids_number=>self%ib%solids_number)
   msg = msg_//' dq '
   call allocate_variable(var=self%dq,       ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%dq = 0._R8P
   msg = msg_//' flx '
   call allocate_variable(var=self%flx,      ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flx = 0._R8P
   msg = msg_//' fly '
   call allocate_variable(var=self%fly,      ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%fly = 0._R8P
   msg = msg_//' flz '
   call allocate_variable(var=self%flz,      ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
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
   call self%adam%io%save_xh5f(basename=trim(output_basename_),    &
                               q=self%field%q, q_name=self%q_name, &
                               with_ghost=with_ghost,              &
                               with_cell_morton=.true.,            &
                               t=self%time%it, time=self%time%time)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_xh5f

   subroutine save_residuals(self)  !invariato ma commentato comando MPI ALLREDUCE
   !< Save residuals history.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call self%field%compute_normL2_residuals(dq=self%dq, norm=self%field%residuals)
      do v=1, self%nv
         !call MPI_ALLREDUCE(MPI_IN_PLACE, self%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
         self%field%residuals(v) = sqrt(self%field%residuals(v))/sqrt(real(self%ni*self%nj*self%nk, R8P))
      enddo
      if (self%mpih%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                           blocks_number=self%blocks_number, residuals=self%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_restart_files(self) !invariato
   !< Save restart files.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self) !ok, commentato parte relativa a q_aux
   !< Save all simulation data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q=self%q)
      !call self%compute_q_auxiliary(q=self%field%q, q_aux=self%q_aux)

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f
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
         elseif (self%physics%D_divergence_cleaner .and. .not.self%physics%B_Divergence_cleaner) then
            call self%slices%save_mat(basename=self%io%output_basename, &
                                      it=self%time%it,                  &
                                      it_max=self%time%it_max,          &
                                      time=self%time%time,              &
                                      time_max=self%time%time_max,      &
                                      adam=self%adam,                   &
                                      q=self%q,                         &
                                      q_name=['Dx  ','Dy  ','Dz  ','Bx  ','By  ','Bz  ','Jx  ','Jy  ','Jz  ','Phi '])
         elseif (self%physics%D_divergence_cleaner .and. self%physics%B_Divergence_cleaner) then
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
   subroutine set_boundary_conditions(self, q) !ok, commentato parte relativa a inflow lasciando estrapolazione e riscritto associate senza R e cv
   !< Set boundary conditions of equation.
   class(prism_cpu_object), intent(in)    :: self              !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:) !< Conservative variables.
   integer(I4P)                        :: b, c, i, j, k, v     !< Counter.
   integer(I4P)                        :: idelta,jdelta,kdelta !< IJK delta step for extrapolation.
   integer(I4P)                        :: bc_type              !< Boundary condition type.
   integer(I4P)                        :: crown                !< Crown counter.
   integer(I4P)                        :: fec                  !< Boundary fec (1 to 26).
   integer(I4P)                        :: fec_1_6              !< Boundary fec (1 to 6).
   integer(I4P)                        :: alfa_D, beta_D, gamma_D !< Indici alfa beta gamma come in Barbas.
   integer(I4P)                        :: alfa_B, beta_B, gamma_B !< Indici alfa beta gamma come in Barbas.
   real(R8P)                           :: s1                   !< Coefficiente pari a +-1.
   real(R8P)                           :: ds                   !< Distanza tra le celle in x, y o z.
   real(R8P)                           :: ngc_r, crown_r           !< Numero di gc totale, reale
   real(R8P)                           :: ref(1:9)             !< Vettore di stato di riferimento per assegnazione gc.
   real(R8P)                           :: fi, f               !< Variabili phi e f fWL.
   !associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown, &
   !          nv=>self%nv, ngc=>self%ngc, cv=>self%physics%eos(1)%cv, R=>self%physics%eos(1)%R, q_bc_vars=>self%bc%q)
   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown, q_old=>self%q_old, &
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

                     !print *, 'i valori del vettore di stato nella gc della faccia', fec, 'sono:'
                     !print *, ref
                     !print *, 'Dx = ', q(1,i,j,k,b)
                     !print *, 'Dy = ', q(2,i,j,k,b)
                     !print *, 'Dz = ', q(3,i,j,k,b)
                     !print *, 'Bx = ', q(4,i,j,k,b)
                     !print *, 'By = ', q(5,i,j,k,b)
                     !print *, 'Bz = ', q(6,i,j,k,b)
                     !print *, 'Jx = ', q(7,i,j,k,b)
                     !print *, 'Jy = ', q(8,i,j,k,b)
                     !print *, 'Jz = ', q(9,i,j,k,b)

                  else
                     do v=1, nv
                       q(v,i,j,k,b) = 0.0_R8P
                     enddo
                  endif

               !elseif (bc_type == BC_INFLOW) then
               !    q(1,i,j,k,b) = q_bc_vars(1, fec_1_6)
               !    q(2,i,j,k,b) = q_bc_vars(1, fec_1_6)* q_bc_vars(2, fec_1_6)
               !    q(3,i,j,k,b) = q_bc_vars(1, fec_1_6)* q_bc_vars(3, fec_1_6)
               !    q(4,i,j,k,b) = q_bc_vars(1, fec_1_6)* q_bc_vars(4, fec_1_6)
               !    q(5,i,j,k,b) = q_bc_vars(1, fec_1_6)*                                &
               !                   (cv*q_bc_vars(5, fec_1_6)/(q_bc_vars(1, fec_1_6)*R) + &
               !                   0.5_R8P*(q_bc_vars(2, fec_1_6)**2+q_bc_vars(3, fec_1_6)**2+q_bc_vars(4, fec_1_6)**2))
               endif
            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self) !ok, resta identico
   !< Set initial conditions and coils on field.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field)

   !< Aggiungo setting iniziale delle spire direttamente dentro IC
   call self%coil%set_coils(physics=self%physics, field=self%field)

   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step) !invariato
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
   subroutine compute_dt(self) !modificata introducendo lmin e ciclo for per calcolo lunghezza minima. commentata parte omp
   class(prism_cpu_object), intent(inout) :: self                            !< The equation.
   real(R8P)                              :: umax                            !< Maximum speed of waves propagation (light speed)
   real(R8P)                              :: lmin                            !< Minimal cell lenght
   real(R8P)                              :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                           :: b, i, j, k                      !< Counter.
   lmin = huge(1._R8P)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, blocks_number=>self%blocks_number, &
            dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), &
            D_divergence_cleaner=>self%physics%D_divergence_cleaner, &
            B_divergence_cleaner=>self%physics%B_divergence_cleaner) !, chi=>self%physics%chi, eta=>self%physics%eta)
   ! !$omp parallel do collapse(4) default(firstprivate) shared(dx,dy,dz,q_aux) reduction(max:umax)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               dx_locale = dx(b)*0.5_R8P
               dy_locale = dy(b)*0.5_R8P
               dz_locale = dz(b)*0.5_R8P
               lmin = min(lmin, dx_locale, dy_locale, dz_locale)
               !umax = max(umax, (abs(q_aux(2,i,j,k,b)) + ss)/dx_locale + 2._R8P*mu/(q_aux(1,i,j,k,b))/dx_locale**2 + &
               !                 (abs(q_aux(3,i,j,k,b)) + ss)/dy_locale + 2._R8P*mu/(q_aux(1,i,j,k,b))/dy_locale**2 + &
               !                 (abs(q_aux(4,i,j,k,b)) + ss)/dz_locale + 2._R8P*mu/(q_aux(1,i,j,k,b))/dz_locale**2)
            enddo
         enddo
      enddo
   enddo
   ! !$omp end parallel do

   if (D_divergence_cleaner) then
      umax = max(self%physics%chi*sqrt(1._R8P/(EPS0*MU0)), self%physics%eta*sqrt(1._R8P/(EPS0*MU0)))
   else
      umax = sqrt(1._R8P/(EPS0*MU0))
   endif

   self%time%dt = self%time%CFL*lmin / umax

   endassociate

    !call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endsubroutine compute_dt

   subroutine compute_residuals(self, q, dq) !tolta da associazione parte physics/eos e i q_aux. Commentato eikonal e q_aux. Tolta parte diffusiva e modificati flussi convettivi per correnti
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
   real(R8P)                              :: vmax   !< Maximum speed of waves propagation.
   integer(I4P)                           :: b, i, j, k !< Counter.

   call self%update_ghost(q=q)

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number, &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), CFL=>self%time%CFL,      &
             phi=>self%ib%phi, flx=>self%flx, fly=>self%fly, flz=>self%flz,                                         &
             weno_s=>self%weno%S,                                                                                   &
             weno_a=>self%weno%a, weno_p=>self%weno%p, weno_d=>self%weno%d, ror_number=>self%weno%ror_number,       &
             ror_schemes=>self%weno%ror_schemes, ror_ivar=>self%weno%ror_ivar,                                      &
             ror_threshold=>self%weno%ror_threshold, enable_ror_stats=>self%weno%enable_ror_stats,                  &
             cell_scheme=>self%weno%cell_scheme, ror_stats=>self%weno%ror_stats, weno_zeps=>self%weno%zeps,         &
             solids_number=>self%ib%solids_number, null_xyz=>self%grid%null_xyz, time=>self%time%time,              &
             A=>self%coil%A, freq=>self%coil%f, phase=>self%coil%phase, coil_flag =>self%coil%coil_flag,            &
             d=>self%coil%d, td=>self%coil%td, chi=>self%physics%chi, eta=>self%physics%eta,                        &
             D_divergence_cleaner=>self%physics%D_divergence_cleaner,                                               &
             B_divergence_cleaner=>self%physics%B_divergence_cleaner)
   if (blocks_number > 0) then
      if (.not.null_xyz(1)) then
         call compute_fluxes_convective(dir=1,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,weno_s=weno_S, &
                                        weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_zeps=weno_zeps,q=q,fluxes=flx,    &
                                        chi=chi, d_divergence_cleaner=d_divergence_cleaner,                              &
                                        b_divergence_cleaner=B_divergence_cleaner)
      else
         call assign_omp(blocks_number=blocks_number, ngc=ngc, lhs=flx, rhs=0._R8P)
      endif
      if (.not.null_xyz(2)) then
         call compute_fluxes_convective(dir=2,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,weno_s=weno_S, &
                                        weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_zeps=weno_zeps,q=q,fluxes=fly,    &
                                        chi=chi, D_divergence_cleaner=D_divergence_cleaner,                              &
                                        B_divergence_cleaner=B_divergence_cleaner)
      else
         call assign_omp(blocks_number=blocks_number, ngc=ngc, lhs=fly, rhs=0._R8P)
      endif
      if (.not.null_xyz(3)) then
         call compute_fluxes_convective(dir=3,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,weno_s=weno_S, &
                                        weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_zeps=weno_zeps,q=q,fluxes=flz,    &
                                        chi=chi, D_divergence_cleaner=D_divergence_cleaner,                              &
                                        B_divergence_cleaner=B_divergence_cleaner)
      else
         call assign_omp(blocks_number=blocks_number, ngc=ngc, lhs=flz, rhs=0._R8P)
      endif

      if (solids_number>0) then
         call compute_fluxes_difference(null_xyz=null_xyz,                                                                   &
                                        blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                        dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, phi=phi, dq=dq, q=q, eta=eta,        &
                                        chi=chi, D_divergence_cleaner=D_divergence_cleaner,                                  &
                                        B_divergence_cleaner=B_divergence_cleaner)
      else
         call compute_fluxes_difference(null_xyz=null_xyz,                                                                   &
                                        blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                        dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, dq=dq, q=q, eta=eta,                 &
                                        chi=chi, D_divergence_cleaner=D_divergence_cleaner,                                  &
                                        B_divergence_cleaner=B_divergence_cleaner)
      endif
   endif
   if (d_divergence_cleaner) then
      vmax = chi/sqrt(EPS0*MU0)
   else
      vmax = sqrt(1._R8P/(EPS0*MU0))
   endif
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               self%field_Div(3,i,j,k,b) = self%coil%coil_flag(i,j,k,b)!CFL/vmax*(dq(1,i,j,k,b)+dq(2,i,j,k,b)+dq(3,i,j,k,b))
               self%field_Div(4,i,j,k,b) = CFL/vmax*(dq(4,i,j,k,b)+dq(5,i,j,k,b)+dq(6,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_residuals

   subroutine integrate(self, do_ghost_syncro) !occhio a come usi il dx, lo hai fatto nel modo più semplice possibile
                                               !ma per griglia non uniforme va modificato
   !< Perform one step integration.
   class(prism_cpu_object), intent(inout)         :: self             !< The equation.
   logical,                 intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   logical                                        :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s,b,i,j,k,var    !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number,   &
             time=>self%time%time, A=>self%coil%A, freq=>self%coil%f, phase=>self%coil%phase,           &
             coil_flag =>self%coil%coil_flag, d=>self%coil%d, td=>self%coil%td, J_vec=>self%coil%J_vec, &
             dx=>self%field%dxyz(1,1))

   if (self%coil%total_coils_number >= 1_I4P) then
      call compute_coils_current(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=self%q, time=time, A=A, d=d, &
                                 f=freq, phase=phase, coil_flag=coil_flag, td=td, J_vec=J_vec, dx=dx)
   endif

   call self%rk%initialize_stages(q=self%q)
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      ! low storage RK working on q_rk_gpu(:,:,:,:,:,1)/q_gpu as stages, update q_gpu in place
      do s=1, self%rk%nrk

         call self%compute_residuals(q=self%q, dq=self%dq)

         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk%compute_stage_ls(s=s,dt=self%time%dt,phi=self%ib%phi,dq=self%dq,q=self%q)
         else
            call self%rk%compute_stage_ls(s=s,dt=self%time%dt,dq=self%dq,q=self%q)
         endif

      enddo
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      ! RK working on q_rk_gpu as stages
      do s=1, self%rk%nrk
         if (self%ib%solids_number>0) then
            call self%rk%compute_stage(s=s, dt=self%time%dt, phi=self%ib%phi)
         else
            call self%rk%compute_stage(s=s, dt=self%time%dt)
         endif
         call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq)
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

   !calcolo della divergenza tramite differenze finite
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
                  !if (self%coil%coil_flag(i,j,k,b) > 0_I4P) then
                  !   self%q(1:6,i,j,k,b) = 0._R8P ! azzero i campi dentro le spire
                  !endif
               self%field_Div(1,i,j,k,b) = 0.5_R8P*((self%q(1,i+1,j,k,b) - self%q(1,i-1,j,k,b))/dx + &
                                          (self%q(2,i,j+1,k,b) - self%q(2,i,j-1,k,b))/dx + &
                                          (self%q(3,i,j,k+1,b) - self%q(3,i,j,k-1,b))/dx)
               self%field_Div(2,i,j,k,b) = 0.5_R8P*((self%q(4,i+1,j,k,b) - self%q(4,i-1,j,k,b))/dx + &
                                          (self%q(5,i,j+1,k,b) - self%q(5,i,j-1,k,b))/dx + &
                                          (self%q(6,i,j,k+1,b) - self%q(6,i,j,k-1,b))/dx)
            enddo
         enddo
      enddo
   enddo

   endassociate
   endsubroutine integrate

   subroutine simulate(self, filename) !invariata ma ho aggiunto parte set coils insieme a ic
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
      !print *, self%q(:,self%ni,1,1,1), 'Riga 927 CPU object'
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

      !self%time%time = self%time%time + self%time%dt
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

   ! non TBP
   subroutine assign_omp_R8P_5D(blocks_number, ngc, lhs, rhs)
   !< Assign array to scalar value with OpenMP threads (kind R8P, rank 5).
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   real(R8P),    intent(inout) :: lhs(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Lest hand side.
   real(R8P),    intent(in)    :: rhs                             !< Right hand side.
   integer(I4P)                :: ni, nj, nk, nv, b, i, j, k, v   !< Counter.

   nv = ubound(lhs,dim=1)
   ni = ubound(lhs,dim=2) - ngc
   nj = ubound(lhs,dim=3) - ngc
   nk = ubound(lhs,dim=4) - ngc
   !!$omp parallel do collapse(5) default(firstprivate) shared(lhs)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do v=1    , nv
      lhs(v,i,j,k,b) = rhs
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine assign_omp_R8P_5D

   subroutine compute_fluxes_convective(dir,blocks_number,ni,nj,nk,ngc,nv,weno_s,weno_a,weno_p,weno_d,weno_zeps,q, &
                                        fluxes,chi,D_divergence_cleaner,B_divergence_cleaner)
   !< Compute convective fluxes along direction `dir`.
   integer(I4P), intent(in)    :: dir                                !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                 !< Number of conservative varibales.
   integer(I4P), intent(in)    :: weno_s                             !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                   !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                          !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)    :: chi                                !< Coefficient to compute transport velocity of field divergence error
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   logical,      intent(in)    :: D_divergence_cleaner               !< Flag to perform electric field divergence cleaning.
   logical,      intent(in)    :: B_divergence_cleaner               !< Flag to perform magnetic field divergence cleaning.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   real(R8P)                   :: el(nv,nv), er(nv,nv)               !< Left and right eigenvalues.
   real(R8P)                   :: fmp (1:2,                   1:nv)  !< Fluxes -+ decomposition.
   real(R8P)                   :: fmpc(1:2,1-weno_s:-1+weno_s,1:nv)  !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:nv)                     !< Fluxes +- reconstructed.
   logical                     :: ror_recompute                      !< Flag to perform ROR.
   integer(I4P)                :: r, v, vv                           !< Counter.
   integer(I4P)                :: b, i, j, k                         !< Counter.
   integer(I4P)                :: si(3), si_i, si_j, si_k            !< Directional (1=x,2=y,3=z) increment.
   real(R8P)                   :: sir(3)                             !< Directional (1=x,2=y,3=z) increment, real.
   !integer(I4P)                :: uni, ut1, ut2                      !< Index of normal and tangential velocities.
   real(R8P)                   :: evmax                              !< Maximum waves speed estimation.
   integer(I4P)                :: s, is, js, ks                      !< Counter.

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

   if (D_divergence_cleaner) then
      evmax = chi*sqrt(1._R8P/(EPS0*MU0))  !velocità massima, sempre pari a quella della luce per i campi ma non per i correttori!
   else
      evmax = sqrt(1._R8P/(EPS0*MU0))
   endif

   !uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   !ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   !ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)

   ! !$omp parallel do collapse(4) default(firstprivate) shared(weno_a, weno_p, weno_d, q_aux, fluxes)
   do b=1, blocks_number
      do k=si_k, nk
         do j=si_j, nj
            do i=si_i, ni
               !call compute_max_eigenvalues(si=si,sir=sir,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,q_aux=q_aux,evmax=evmax)
               call decompose_fluxes_convective_llf(sir=sir, nv=nv, q=q(:,i      ,j      ,k      ,b), evmax=evmax, chi=chi, &
                                                    fmp=fmp, B_divergence_cleaner=B_divergence_cleaner, &
                                                    D_divergence_cleaner=D_divergence_cleaner)
               call decompose_fluxes_convective_llf(sir=sir, nv=nv, q=q(:,i+si(1),j+si(2),k+si(3),b), evmax=evmax, chi=chi, &
                                                    fmp=fpmr, B_divergence_cleaner=B_divergence_cleaner, &
                                                    D_divergence_cleaner=D_divergence_cleaner)
               fluxes(:,i,j,k,b) = fmp(2,:) + fpmr(1,:)
            enddo
         enddo
      enddo
   enddo
   !!$omp end parallel do
   endsubroutine compute_fluxes_convective

   subroutine compute_fluxes_convective_mod(dir,blocks_number,ni,nj,nk,ngc,nv,weno_s,weno_a,weno_p,weno_d,weno_zeps,q, &
               fluxes,chi,D_divergence_cleaner,B_divergence_cleaner, coil_flag, dx, dy, dz) !cambiato q_aux con q, tolto g dagli input, commentata openmp
   !< Compute convective fluxes along direction `dir`.
   !< Versione modificata per approssimare la J come corrente superficiale e vedere se va meglio l'errore du D e il raggiungimento
   !< dello stato stazionario.
   integer(I4P), intent(in)    :: dir                                         !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                               !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                          !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                          !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                          !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                         !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                          !< Number of conservative varibales.
   integer(I4P), intent(in)    :: weno_s                                      !< Weno stencils number/dimension.
   integer(I4P), intent(in)    :: coil_flag(1-ngc:,1-ngc:,1-ngc:,1:)      !< Matrice contenente informazioni su quale spira pass per una certa cella
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                            !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                         !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                         !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                                   !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)    :: chi                                         !< Coefficient to compute transport velocity of field divergence error
   real(R8P),    intent(in)    :: dx(1:), dy(1:), dz(1:)                      !< Grid cell sizes.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)               !< Field variables.
   logical,      intent(in)    :: D_divergence_cleaner                        !< Flag to perform electric field divergence cleaning.
   logical,      intent(in)    :: B_divergence_cleaner                        !< Flag to perform magnetic field divergence cleaning.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:)          !< Fluxes.
   real(R8P)                   :: el(nv,nv), er(nv,nv)                        !< Left and right eigenvalues.
   real(R8P)                   :: fmp (1:2,                   1:nv)           !< Fluxes -+ decomposition.
   real(R8P)                   :: fmpc(1:2,1-weno_s:-1+weno_s,1:nv)           !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:nv)                              !< Fluxes +- reconstructed.
   logical                     :: ror_recompute                               !< Flag to perform ROR.
   integer(I4P)                :: r, v, vv                                    !< Counter.
   integer(I4P)                :: b, i, j, k                                  !< Counter.
   integer(I4P)                :: si(3), si_i, si_j, si_k                     !< Directional (1=x,2=y,3=z) increment.
   real(R8P)                   :: sir(3)                                      !< Directional (1=x,2=y,3=z) increment, real.
   real(R8P)                   :: vecR(1:nv),vecL(1:nv)                       !< Vettori di appoggio per gli stati del problema di riemann.
   !integer(I4P)                :: uni, ut1, ut2                              !< Index of normal and tangential velocities.
   real(R8P)                   :: evmax                                       !< Maximum waves speed estimation.
   integer(I4P)                :: s, is, js, ks                               !< Counter.

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

   if (D_divergence_cleaner) then
      evmax = chi*sqrt(1._R8P/(EPS0*MU0))  !velocità massima, sempre pari a quella della luce per i campi ma non per i correttori!
   else
      evmax = sqrt(1._R8P/(EPS0*MU0))
   endif

   !uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   !ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   !ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)

   ! !$omp parallel do collapse(4) default(firstprivate) shared(weno_a, weno_p, weno_d, q_aux, fluxes)
   do b=1, blocks_number
      do k=si_k, nk
         do j=si_j, nj
            do i=si_i, ni
               !call compute_max_eigenvalues(si=si,sir=sir,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,q_aux=q_aux,evmax=evmax)
               if (coil_flag(i,j,k,b) > 0_I4P) then
                  vecL(1) = -q(1,i+si(1),j+si(2),k+si(3),b)
                  vecL(2) = -q(2,i+si(1),j+si(2),k+si(3),b)
                  vecL(3) = -q(3,i+si(1),j+si(2),k+si(3),b)
                  vecL(4:6) = q(4:6,i+si(1),j+si(2),k+si(3),b)-MU0*crossproduct(sir,q(7:9,i,j,k,b))*dx(b)
                  vecL(7:9) = q(7:9,i,j,k,b)
                  !print *, vecL, 'Riga 1039 CPU object'
               else
                  vecL = q(:,i,j,k,b)
               endif
               call decompose_fluxes_convective_llf(sir=sir, nv=nv, q=vecL, evmax=evmax, chi=chi, fmp=fmp, &
                                                    B_divergence_cleaner=B_divergence_cleaner, &
                                                    D_divergence_cleaner=D_divergence_cleaner)
               if (coil_flag(i+si(1),j+si(2),k+si(3),b) > 0_I4P) then
                  vecR(1) = -q(1,i,j,k,b)
                  vecR(2) = -q(2,i,j,k,b)
                  vecR(3) = -q(3,i,j,k,b)
                  vecR(4:6) = q(4:6,i,j,k,b)+MU0*crossproduct(sir,q(7:9,i+si(1),j+si(2),k+si(3),b))*dx(b)
                  vecR(7:9) = q(7:9,i+si(1),j+si(2),k+si(3),b)
                  !print *, vecR, 'Riga 1040 CPU object'
               else
                  vecR = q(:,i+si(1),j+si(2),k+si(3),b)
               endif

               call decompose_fluxes_convective_llf(sir=sir, nv=nv, q=vecR, evmax=evmax, chi=chi, fmp=fpmr, &
                                                    B_divergence_cleaner=B_divergence_cleaner, &
                                                    D_divergence_cleaner=D_divergence_cleaner)
               fluxes(:,i,j,k,b) = fmp(2,:) + fpmr(1,:)
               if (coil_flag(i,j,k,b) > 0_I4P .and. coil_flag(i+si(1),j+si(2),k+si(3),b) > 0_I4P) then
                  fluxes(:,i,j,k,b) = 0._R8P
               endif
            enddo
         enddo
      enddo
   enddo
   !!$omp end parallel do
   endsubroutine compute_fluxes_convective_mod

   subroutine compute_fluxes_difference(null_xyz, blocks_number, ni, nj, nk, ngc, nv, ib_eps, dx, dy, dz, flx, fly, flz, &
                                       phi, dq, q, eta, chi, D_divergence_cleaner, B_divergence_cleaner)
                                       !commentata parte openmp e aggiunto vettore di stato agli input per sfruttare correnti
                                       !< Compute fluxes difference.
   logical,      intent(in)           :: null_xyz(3)                     !< Nullified directions tags.
   logical,      intent(in)           :: D_divergence_cleaner            !< Flag to perform electric field divergence cleaning.
   logical,      intent(in)           :: B_divergence_cleaner            !< Flag to perform magnetic field divergence cleaning.
   integer(I4P), intent(in)           :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)           :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                              !< Number of conservative varibales.
   real(R8P),    intent(in)           :: chi, eta                        !< Coefficiente modello correzione divergenza campi.
   real(R8P),    intent(in)           :: ib_eps                          !< Tolerance IB delta ratio.
   real(R8P),    intent(in)           :: dx(1:), dy(1:), dz(1:)          !< Space steps.
   real(R8P),    intent(in)           :: flx(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)           :: fly(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)           :: flz(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in)           :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< State variables vector
   real(R8P),    intent(in), optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),    intent(inout)        :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes differences.
   real(R8P)                          :: delta_x, delta_y, delta_z       !< Space steps.
   real(R8P)                          :: dx_locale, dy_locale, dz_locale !< Local space steps.
   real(R8P)                          :: qmx, qmy, qmz                   !< Momentum nullification scalar.
   integer(I4P)                       :: b, i, j, k, v                   !< Counter.
   integer(I4P)                       :: all_solids                      !< Last phi index, all solids summary.

   qmx = 1._R8P ; if (null_xyz(1)) qmx = 0._R8P
   qmy = 1._R8P ; if (null_xyz(2)) qmy = 0._R8P
   qmz = 1._R8P ; if (null_xyz(3)) qmz = 0._R8P
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
   !    !$omp parallel do collapse(4) default(firstprivate) shared(dx,dy,dz,flx,fly,flz,phi,dq)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         dx_locale = dx(b)
         if (phi(all_solids,i,j,k,b)<0.) then
            if (phi(all_solids,i+1,j,k,b)*phi(all_solids,i-1,j,k,b)<0) then
               if (phi(all_solids,i+1,j,k,b)>0.) then
                  delta_x = -phi(all_solids,i,j,k,b)/(phi(all_solids,i+1,j,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dx(b)
                  dx_locale = dx(b)/2 + delta_x
               else
                  delta_x = -phi(all_solids,i,j,k,b)/(phi(all_solids,i-1,j,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dx(b)
                  dx_locale = dx(b)/2 + delta_x
               endif
            endif
         endif
         dy_locale = dy(b)
         if (phi(all_solids,i,j,k,b)<0.) then
            if (phi(all_solids,i,j+1,k,b)*phi(all_solids,i,j-1,k,b)<0) then
               if (phi(all_solids,i,j+1,k,b)>0.) then
                  delta_y = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j+1,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dy(b)
                  dy_locale = dy(b)/2 + delta_y
               else
                  delta_y = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j-1,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dy(b)
                  dy_locale = dy(b)/2 + delta_y
               endif
            endif
         endif
         dz_locale = dz(b)
         if (phi(all_solids,i,j,k,b)<0.) then
            if (phi(all_solids,i,j,k+1,b)*phi(all_solids,i,j,k-1,b)<0) then
               if (phi(all_solids,i,j,k+1,b)>0.) then
                  delta_z = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j,k+1,b)-phi(all_solids,i,j,k,b)+ib_eps)*dz(b)
                  dz_locale = dz(b)/2 + delta_z
               else
                  delta_z = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j,k-1,b)-phi(all_solids,i,j,k,b)+ib_eps)*dz(b)
                  dz_locale = dz(b)/2 + delta_z
               endif
            endif
         endif
         do v=1, nv
            dq(v,i,j,k,b) = - (flx(v,i,j,k,b)-flx(v,i-1,j,k,b))/dx_locale &
                            - (fly(v,i,j,k,b)-fly(v,i,j-1,k,b))/dy_locale &
                            - (flz(v,i,j,k,b)-flz(v,i,j,k-1,b))/dz_locale
         enddo
         !dq(2,i,j,k,b) = dq(2,i,j,k,b) * qmx
         !dq(3,i,j,k,b) = dq(3,i,j,k,b) * qmy
         !dq(4,i,j,k,b) = dq(4,i,j,k,b) * qmz

         !Completo calcolo aggiungendo termini sorgenti legato alle correnti delle spire (per ora)
         !E ad una eventuale correzione parabolica della divergenza (parametro eta diverso da zero)

         dq(1,i,j,k,b) = dq(1,i,j,k,b) - q(7,i,j,k,b)
         dq(2,i,j,k,b) = dq(2,i,j,k,b) - q(8,i,j,k,b)
         dq(3,i,j,k,b) = dq(3,i,j,k,b) - q(9,i,j,k,b)

         if (D_divergence_cleaner .and. .not.B_divergence_cleaner .and. eta>0._R8P) then
            dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
         elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
            dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
            dq(11,i,j,k,b) = dq(11,i,j,k,b) - chi/eta*chi/eta*q(11,i,j,k,b)
         endif

      enddo
      enddo
      enddo
      enddo
   !    !$omp end parallel do
   else
   !    !$omp parallel do collapse(4) default(firstprivate) shared(dx,dy,dz,flx,fly,flz,phi,dq)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         do v=1, nv
            dq(v,i,j,k,b) = - (flx(v,i,j,k,b)-flx(v,i-1,j,k,b))/dx(b) &
                            - (fly(v,i,j,k,b)-fly(v,i,j-1,k,b))/dy(b) &
                            - (flz(v,i,j,k,b)-flz(v,i,j,k-1,b))/dz(b)
         enddo
         !dq(2,i,j,k,b) = dq(2,i,j,k,b) * qmx
         !dq(3,i,j,k,b) = dq(3,i,j,k,b) * qmy
         !dq(4,i,j,k,b) = dq(4,i,j,k,b) * qmz

         !Completo calcolo aggiungendo termini sorgenti legato alle correnti delle spire (per ora)
         !E ad una eventuale correzione parabolica della divergenza (parametro eta diverso da zero)

         dq(1,i,j,k,b) = dq(1,i,j,k,b) - q(7,i,j,k,b)
         dq(2,i,j,k,b) = dq(2,i,j,k,b) - q(8,i,j,k,b)
         dq(3,i,j,k,b) = dq(3,i,j,k,b) - q(9,i,j,k,b)

         if (D_divergence_cleaner .and. .not.B_divergence_cleaner .and. eta>0._R8P) then
            dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
         elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
            dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
            dq(11,i,j,k,b) = dq(11,i,j,k,b) - chi/eta*chi/eta*q(11,i,j,k,b)
         endif

      enddo
      enddo
      enddo
      enddo
   !    !$omp end parallel do
   endif
   endsubroutine compute_fluxes_difference

   subroutine decompose_fluxes_convective_llf(sir, nv, q, evmax, chi, fmp, B_divergence_cleaner, D_divergence_cleaner)
   !< Decompose convective fluxes using the Local-Lax-Friedrichs (LLF, Rusanov) approximation
   logical,      intent(in)    :: D_divergence_cleaner   !< Flag to perform electric field divergence cleaning.
   logical,      intent(in)    :: B_divergence_cleaner   !< Flag to perform magnetic field divergence cleaning.
   integer(I4P), intent(in)    :: nv                     !< Number of conservative varibales.
   real(R8P),    intent(in)    :: sir(3)                 !< Directional (1=x,2=y,3=z) increment, real
   !real(R8P),    intent(in)    :: q_aux(1:)             !< Auxiliary variables.
   real(R8P),    intent(in)    :: evmax                  !< Maximum waves speeds estimation.
   real(R8P),    intent(in)    :: chi                    !< Coefficient to compute transport velocity of field divergence error
   real(R8P),    intent(in)    :: q(1:)                  !< Conservative variables.
   real(R8P),    intent(inout) :: fmp(1:,1:)             !< Fluxes, negative/positive terms [1:2,1:nv].
   real(R8P)                   :: f(1:nv)                !< Conservative fluxes.
   integer(I4P)                :: v                      !< Counter.

   !call compute_conservatives_scalar(q_aux=q_aux,q=q)
   !call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux,f=f)
   if (.not.D_divergence_cleaner .and. .not.B_divergence_cleaner) then
      call compute_convective_fluxes_Maxwell(sir=sir,q=q,f=f)
   elseif (D_divergence_cleaner .and. .not.B_divergence_cleaner) then
      call compute_convective_fluxes_Maxwell_div_D(sir=sir,q=q,f=f,chi=chi)
   elseif (D_divergence_cleaner .and. B_divergence_cleaner) then
      call compute_convective_fluxes_Maxwell_div_D_B(sir=sir,q=q,f=f,chi=chi)
   endif

   do v=1, nv !(nv-3_I4P)
      fmp(2,v) = 0.5_R8P * (f(v) + evmax * q(v))
      fmp(1,v) = f(v) - fmp(2,v)
   enddo
      fmp(:,7) = 0._R8P
      fmp(:,8) = 0._R8P
      fmp(:,9) = 0._R8P
   endsubroutine decompose_fluxes_convective_llf

   subroutine compute_coils_current(ni, nj, nk, ngc, blocks_number, q, time, A, d, f, phase, coil_flag, td, J_vec, dx)

      integer(I4P), intent(in)           :: blocks_number                               !< Number of blocks.
      integer(I4P), intent(in)           :: ni                                          !< Grid cells number in I direction.
      integer(I4P), intent(in)           :: nj                                          !< Grid cells number in J direction.
      integer(I4P), intent(in)           :: nk                                          !< Grid cells number in K direction.
      integer(I4P), intent(in)           :: ngc                                         !< Ghost cells number.
      integer(I4P), intent(in)           :: coil_flag(1-ngc:,1-ngc:,1-ngc:,1:)      !< Matrice contenente informazioni su quale spira pass per una certa cella
      real(R8P),    intent(in)           :: time                                        !< Simulation time, to compute current value if AC
      real(R8P),    intent(in)           :: dx                                          !< Space step in x direction (m)
      real(R8P),    intent(in)           :: A(1:)                                       !< Current amplitude (A)
      real(R8P),    intent(in)           :: f(1:)                                       !< Current frequency, if AC (Hz)
      real(R8P),    intent(in)           :: phase(1:)                                   !< Current initial phase, if AC
      real(R8P),    intent(in)           :: d(1:)                                       !< Wire diameter
      real(R8P),    intent(in)           :: td                                          !< Delay di accensione della spira
      real(R8P),    intent(in)           :: J_vec(1:,1-ngc:,1-ngc:,1-ngc:,1:)           !< Matrice versori di corrente delle spire nelle celle
      real(R8P),    intent(inout)        :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)               !< Field variables.
      real(R8P)                          :: current_density                             !< Current density
      real(R8P)                          :: g                                           !< Polinomio caratteristico transitorio accensione spira
      integer(I4P)                       :: coil_id                                     !< ID per identificare spira
      integer(I4P)                       :: i,j,k,b,n                                   !< Counter

      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  coil_id = coil_flag(i,j,k,b)
                  if (coil_id > 0_I4P) then
                     !Per DC frequenza e fase sono nulle, quindi se uso la funzione coseno
                     !mi rispramio anche il selectcase

                     !Densità di corrente al tempo t della spira n-esima identificata da (coil_id)
                     !current_density = 4*A(coil_id)/(pi*d(coil_id)**2)*cos(2*pi*f(coil_id)*time + phase(coil_id)*pi/180.0_R8P)

                     !Modifico calcolo densità di corrente considerando sezione quadrata, per coerenza con calcolo Filippo
                     !E aggiungo transitorio di corrente
                     if (time < td) then
                        g = 10._R8P*(time/td)**3 - 15._R8P*(time/td)**4 + 6._R8P*(time/td)**5
                        current_density = g*A(coil_id)/((d(coil_id)-dx)**2)*cos(phase(coil_id)*pi/180.0_R8P)
                     else
                        current_density = A(coil_id)/((d(coil_id)-dx)**2)*cos(2*pi*f(coil_id)*(time-td) + &
                        phase(coil_id)*pi/180.0_R8P)
                     endif

                     q(7,i,j,k,b) = current_density* J_vec(1,i,j,k,b)
                     q(8,i,j,k,b) = current_density* J_vec(2,i,j,k,b)
                     q(9,i,j,k,b) = current_density* J_vec(3,i,j,k,b)

                     !if (sq_norm(q(7:9,i,j,k,b)) == 0._R8P) then
                        !Se la densità di corrente è nulla non faccio nulla
                     !   q(7:9,i,j,k,b) = current_density*q(7:9,i,j,k,b)
                     !else
                        !Se la densità di corrente è diversa da zero allora rinormalizzo il vettore corrente
                        !e lo moltiplico per la densità di corrente
                     !   q(7,i,j,k,b) = q(7,i,j,k,b)/(sq_norm(q(7:9,i,j,k,b)))**0.5 !lo devo rinormalizzare ogni volta
                     !   q(8,i,j,k,b) = q(8,i,j,k,b)/(sq_norm(q(7:9,i,j,k,b)))**0.5
                     !   q(9,i,j,k,b) = q(9,i,j,k,b)/(sq_norm(q(7:9,i,j,k,b)))**0.5
                     !   q(7:9,i,j,k,b) = current_density*q(7:9,i,j,k,b)
                     !endif
                  endif
               enddo
            enddo
         enddo
      enddo
   endsubroutine compute_coils_current

   function sq_norm(a) result(sq)
   !< Return the square of the norm of vector.
   !<
   !< The square norm if defined as \( N = x^2  + y^2  + z^2 \).
   !<
   !<```fortran
   !< type(vector_RPP) :: pt
   !< pt = ex_RPP + ey_RPP
   !< print "(F3.1)", pt%sq_norm()
   !<```
   !=> 2.0 <<<
   !<
   !<```fortran
   !< type(vector_RPP) :: pt
   !< pt = ex_RPP + ey_RPP
   !< print "(F3.1)", sq_norm_RPP(pt)
   !<```
   !=> 2.0 <<<
   real(R8P), intent(in)  :: a(3)     !< Input vector
   real(R8P)              :: sq       !< Square norm of input

   sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
   endfunction sq_norm

   function crossproduct(a, b) result(cross)

   !< Compute the cross product.
   !<
   !< $$ \vec V=\left({y_1 z_2 - z_1 y_2}\right)\vec i +
   !<           \left({z_1 x_2 - x_1 z_2}\right)\vec j +
   !<           \left({x_1 y_2 - y_1 x_2}\right)\vec k $$
   !< where \( x_i \), \( y_i \) and \( z_i \) \( i=1,2 \) are the components of the vectors.
   !<
   !<```fortran
   !< type(vector_RPP) :: pt(0:2)
   !< pt(1) = 2 * ex_RPP
   !< pt(2) = ex_RPP
   !< pt(0) = pt(1).cross.pt(2)
   !< print "(3(F3.1,1X))", abs(pt(0)%x), abs(pt(0)%y), abs(pt(0)%z)
   !<```
   !=> 0.0 0.0 0.0 <<<

   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))

   endfunction crossproduct
endmodule adam_prism_cpu_object

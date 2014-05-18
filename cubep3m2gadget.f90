program test
  implicit none
  integer(4) :: i, np_local, nts, cur_checkpoint, cur_projection, cur_halofind
  real(4) :: a, t, tau, dt_f_acc, dt_pp_acc, dt_c_acc, mass_p
  real(4), allocatable, dimension(:,:) :: xv
  integer(8), allocatable, dimension(:) :: PID
  integer(8) :: mpi_npart(6), mpi_nparttotal(6)
  integer(4) :: g_npart(6),g_flag_sfr,g_flag_feedback,g_npartTotal(6),g_flag_cooling,g_num_files
  integer(4) :: g_flag_stellarage, g_flag_metals, g_nhighword(6), g_filler(16)
  real(8) :: g_mass(6), g_time, g_redshift, g_Boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam
  real(4) :: boxsize, Omega0, OmegaLambda, HubbleParam
  integer(4) :: num_files
#define EXTRAPID

#ifdef EXTRAPID
  open(unit=21,file="/scratch/00506/ilievit/cubepm_130315_6_1728_47Mpc_ext2/results/8.064xv0.dat",status='old',form='binary')
  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
       cur_projection,cur_halofind,mass_p

  allocate(PID(1:np_local))
  allocate(xv(1:6,1:np_local))

  read(21) xv(:,:np_local)
  close(21)

  open(unit=21,file="/scratch/00506/ilievit/cubepm_130315_6_1728_47Mpc_ext2/results/8.064PID0.dat",status='old',form='binary')
  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
       cur_projection,cur_halofind,mass_p

  read(21) PID
  close(21)
#endif

  OmegaLambda = 0.7
  Omega0 = 0.3
  HubbleParam = 0.7
  num_files = 1  

  g_mass(1:6) = 0.0d0
  g_mass(2) = real(mass_p,8)
  g_npart(1:6) = 0
  g_npart(2) = np_local
  g_npartTotal = g_npart
  g_flag_sfr = 0
  g_flag_feedback = 0
  g_flag_cooling = 0
  g_num_files = num_files
  g_Boxsize = boxsize
  g_Omega0 = Omega0
  g_OmegaLambda = OmegaLambda
  g_HubbleParam = HubbleParam

  open(unit=21,file="../test.bin.0",form='unformatted')
  
  write(21) g_npart, g_mass, g_time, g_redshift, g_flag_sfr, g_flag_feedback, g_npartTotal, &
       g_flag_cooling, g_num_files, g_boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam, &
       g_flag_stellarage, g_flag_metals, g_nhighword, g_filler
  write(21) xv(1:3,1:np_local)
  write(21) xv(4:6,1:np_local)
  write(21) PID(1:np_local)
  close(21)
  
  deallocate(xv,PID)
end program test

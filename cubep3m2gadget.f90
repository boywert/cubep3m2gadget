program test
  use mpi
  implicit none
  integer(4) :: highword,lowword
  integer(4) :: np_local, nts, cur_checkpoint, cur_projection, cur_halofind
  real(4) :: a, t, tau, dt_f_acc, dt_pp_acc, dt_c_acc, mass_p
  real(4), allocatable, dimension(:,:) :: xv
  integer(8), allocatable, dimension(:) :: PID
  integer(8) :: mpi_npart, mpi_nparttotal
  integer(4) :: g_npart(6),g_flag_sfr,g_flag_feedback,g_npartTotal(6),g_flag_cooling,g_num_files
  integer(4) :: g_flag_stellarage, g_flag_metals, g_nhighword(6), g_filler(16)
  real(8) :: g_mass(6), g_time, g_redshift, g_Boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam
  real(4) :: boxsize, Omega0, OmegaLambda, HubbleParam
  integer(4) :: num_files
  character(len=100) :: str_rank,z_s,xv_input,pid_input,output
  real(8) :: redshift
  integer(4) :: totalnodes,rank,ierr
  integer(4) :: ngdim,npdim,ncdim
  real(4) :: nc_offset(3)
  integer(4) :: i,j,k
  real(4) :: gadgetmass

  call mpi_init(ierr)
  npdim = 1728
  ngdim = 6 
  OmegaLambda = 0.73
  Omega0 = 0.27
  HubbleParam = 0.7
  num_files = 216
  boxsize = 47.0         !Mpc/h
  redshift = 8.064
  gadgetmass = 1.0e10    !Msun


  ncdim = 2*npdim
  call mpi_comm_size(mpi_comm_world,totalnodes,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
  ! if(totalnodes /= num_files) then
  !    if(rank == 0) then
  !       print*, "total nodes:",totalnodes,"total files:",num_files
  !       print*, "aborting"
  !       call abort
  !    endif
  ! endif
  write(z_s, "(f10.3)") redshift
  z_s = adjustl(z_s)
  write(str_rank, "(I10)") rank
  str_rank = adjustl(str_rank)
  xv_input = "/scratch/00506/ilievit/cubepm_130315_6_1728_47Mpc_ext2/results/"//z_s(1:len_trim(z_s))//"xv"//str_rank(1:len_trim(str_rank))//".dat"
  pid_input = "/scratch/00506/ilievit/cubepm_130315_6_1728_47Mpc_ext2/results/"//z_s(1:len_trim(z_s))//"PID"//str_rank(1:len_trim(str_rank))//".dat"
  output = "/scratch/01937/cs390/test/"//z_s(1:len_trim(z_s))//"/"
  call system("mkdir -p "//trim(output))
  output = trim(output)//"/"//z_s(1:len_trim(z_s))//"."//str_rank(1:len_trim(str_rank))

  do k=0,ngdim-1
     do j=0,ngdim-1
        do i=0,ngdim-1
           if(rank == i+j*ngdim+k*ngdim**2) then
              nc_offset(1) = real(i*ncdim/ngdim)
              nc_offset(2) = real(j*ncdim/ngdim)
              nc_offset(3) = real(k*ncdim/ngdim)
           endif
        enddo
     enddo
  enddo
#define EXTRAPID

#ifdef EXTRAPID
  open(unit=21,file=trim(xv_input),status='old',form='binary')
  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
       cur_projection,cur_halofind,mass_p

  allocate(PID(1:np_local))
  allocate(xv(1:6,1:np_local))

  read(21) xv(:,:np_local)
  close(21)

  open(unit=21,file=trim(pid_input),status='old',form='binary')
  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
       cur_projection,cur_halofind,mass_p
  if (rank==0) print*, np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
       cur_projection,cur_halofind,mass_p


  read(21) PID
  close(21)
  ! convert cubep3m units -> gadget units 
  do i=1,3
     xv(i,1:np_local) = ((xv(i,1:np_local) + nc_offset(i)))*boxsize/real(ncdim)
  enddo
  print*, xv(4:6,1:np_local)
  !xv(4:6,1:np_local) = xv(4:6,1:np_local)/sqrt(a)

#endif

  g_time = a

  g_mass(1:6) = 0.0d0
  g_mass(2) = real(mass_p,8)
  g_npart(1:6) = 0
  g_npart(2) = np_local
  mpi_npart = int(np_local,8)

  g_flag_sfr = 0
  g_flag_feedback = 0
  g_flag_cooling = 0
  g_num_files = num_files
  g_Boxsize = boxsize
  g_Omega0 = Omega0
  g_redshift = redshift
  g_OmegaLambda = OmegaLambda
  g_HubbleParam = HubbleParam

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call mpi_allreduce(mpi_npart,mpi_nparttotal,1,mpi_integer8,mpi_sum,mpi_comm_world,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  g_npartTotal(1:6) = 0
  g_nhighword(1:6) = 0
  if(rank == 0 ) print*,"total N",mpi_nparttotal
  g_npartTotal(2) = iand(mpi_nparttotal,2**32-1)
  if(rank == 0) print*, g_npartTotal
  g_nhighword(2) = ishft(mpi_nparttotal,-32)
  if(rank == 0) print*,g_nhighword
  open(unit=21,file=trim(output),form='unformatted')

  write(21) g_npart, g_mass, g_time, g_redshift, g_flag_sfr, g_flag_feedback, g_npartTotal, &
       g_flag_cooling, g_num_files, g_boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam, &
       g_flag_stellarage, g_flag_metals, g_nhighword, g_filler
  write(21) xv(1:3,1:np_local)
  write(21) xv(4:6,1:np_local)
  write(21) PID(1:np_local)
  close(21)

  deallocate(xv,PID)

  call mpi_finalize(ierr)
end program test

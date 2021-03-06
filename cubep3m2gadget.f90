module mympi
  use mpi 
  implicit none
  integer(4) :: totalnodes,rank,ierr    
  contains 
    subroutine mympi_init()
      implicit none
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,totalnodes,ierr)
      if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
      call mpi_comm_rank(mpi_comm_world,rank,ierr)
      if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
      if(rank == 0) print*,"Total MPI Tasks:",totalnodes
    end subroutine mympi_init
    subroutine mympi_end()
      implicit none
      call mpi_finalize(ierr)
    end subroutine mympi_end
end module mympi

module cube2gadget_module
  use mympi
  use mpi
  implicit none
  integer(4), parameter :: maxsnap = 1000
  integer(4) :: totalsnaps
  real(4), parameter :: H0 = 100.   ![h*km]/[sec*Mpc]
  real(4), parameter :: RHO_CRIT_0 = 2.7755397e11   ! [h^2*Msun]/[Mpc^3]27.7473406739
  real(4) :: Omega0,OmegaLambda,HubbleParam,boxsize
  real(4) :: c_vunit, c_munit, c_lunit
  integer(4) :: ngdim,npdim,ncdim,num_files
  real(4) :: redshift_list(maxsnap)
  contains
    subroutine cubep3m_config_init(redshift)
      implicit none
      real(4) :: redshift
      num_files = ngdim**3
      ncdim = 2*npdim
      c_vunit = vunit_compute(redshift)
      c_munit = munit_compute()
      c_lunit = lunit_compute()
      return
    end subroutine cubep3m_config_init

    function vunit_compute(redshift)
      real(4) :: vunit_compute
      real(4) :: redshift
      vunit_compute = boxsize * 1.5 * sqrt(Omega0) * H0 / (ncdim ) * (1.+redshift)  !km/s
      return
    end function vunit_compute

    function munit_compute()
      real(4) :: munit_compute
      munit_compute = boxsize * boxsize * boxsize * Omega0 * RHO_CRIT_0 / real(int(npdim,8)**3) / 1.0e10  ! 1.e10 Msun/h
      return 
    end function munit_compute


    function lunit_compute()
      real(4) :: lunit_compute
      lunit_compute = boxsize/real(ncdim) ! Mpc/h
      return 
    end function lunit_compute

    subroutine cube2gadget_call(input_prefix,output_prefix,name_prefix,isnap)
      implicit none
      character(len=1000) :: input_prefix,output_prefix,name_prefix
      integer(4) :: highword,lowword
      integer(4) :: np_local, nts, cur_checkpoint, cur_projection, cur_halofind
      real(4) :: a, t, tau, dt_f_acc, dt_pp_acc, dt_c_acc, mass_p
      real(4), allocatable, dimension(:,:) :: xv
      integer(8), allocatable, dimension(:) :: PID
      integer(8) :: mpi_npart, mpi_nparttotal
      integer(4) :: g_npart(6),g_flag_sfr,g_flag_feedback,g_npartTotal(6),g_flag_cooling,g_num_files
      integer(4) :: g_flag_stellarage, g_flag_metals, g_nhighword(6), g_filler(16)
      real(8) :: g_mass(6), g_time, g_redshift, g_Boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam

      character(len=1000) :: str_rank,z_s,xv_input,pid_input,output,numsnap
      real(4) :: redshift
      real(4) :: nc_offset(3)
      integer(4) :: isnap,iz,i,j,k
      integer(4) :: offset_rank
      iz = isnap - 1
      offset_rank = 0
      redshift = redshift_list(isnap)
      if(rank == 0) print*, "redshift:",redshift
      if(redshift < 0.0) return
      call cubep3m_config_init(redshift)


      ! if(totalnodes /= num_files) then
      !    if(rank == 0) then
      !       print*, "total nodes:",totalnodes,"total files:",num_files
      !       print*, "aborting"
      !       call abort
        !    endif
      ! endif


      write(z_s, "(f10.3)") redshift
      z_s = adjustl(z_s)
      write(str_rank, "(I10)") rank+offset_rank
      str_rank = adjustl(str_rank)
      write(numsnap, "(I3.3)") iz
      numsnap = adjustl(numsnap)
      input_prefix = adjustl(input_prefix)
      output_prefix = adjustl(output_prefix)
#ifndef IC
      xv_input = input_prefix(1:len_trim(input_prefix))//z_s(1:len_trim(z_s))//"xv"//str_rank(1:len_trim(str_rank))//".dat"
      pid_input = input_prefix(1:len_trim(input_prefix))//z_s(1:len_trim(z_s))//"PID"//str_rank(1:len_trim(str_rank))//".dat"
       output = output_prefix(1:len_trim(output_prefix))//"snapdir_"// numsnap(1:len_trim(numsnap))//"/" 
#else
      xv_input = input_prefix(1:len_trim(input_prefix))//z_s(1:len_trim(z_s))//"xv"//str_rank(1:len_trim(str_rank))//".IC"
      pid_input = input_prefix(1:len_trim(input_prefix))//z_s(1:len_trim(z_s))//"PID"//str_rank(1:len_trim(str_rank))//".IC"
      output = output_prefix(1:len_trim(output_prefix))//"ICs/" 
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank == 0) call system("mkdir -p "//trim(output))
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      name_prefix = adjustl(name_prefix)
#ifndef IC
      output = trim(output)//"/"//name_prefix(1:len_trim(name_prefix))//"_"//numsnap(1:len_trim(numsnap))//"."//str_rank(1:len_trim(str_rank))
#else
      output = trim(output)//"/"//name_prefix(1:len_trim(name_prefix))//"_ics.gdt."//str_rank(1:len_trim(str_rank))
#endif


#define EXTRAPID
#ifdef EXTRAPID
      open(unit=21,file=trim(xv_input),status='old',form='binary')
#ifndef IC
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
           cur_projection,cur_halofind,mass_p
#else
      read(21) np_local
#endif
      allocate(PID(1:np_local))
      allocate(xv(1:6,1:np_local))

      read(21) xv(:,:np_local)
      close(21)

      open(unit=21,file=trim(pid_input),status='old',form='binary')
#ifndef IC
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
           cur_projection,cur_halofind,mass_p
      if (rank==0) print*, np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
           cur_projection,cur_halofind,mass_p
#else
      read(21) np_local
      if (rank==0) print*, np_local
#endif
      


      read(21) PID
      close(21)

      if(rank == 0) then
         print*, "lunit:",c_lunit,"Mpc/h"
         print*, "vunit:",c_vunit,"km/s"
         print*, "munit:",c_munit,"1.0e10 Msun/h"
      endif
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
      ! convert cubep3m units -> gadget units 
      do i=1,3
         xv(i,1:np_local) = ((xv(i,1:np_local) + nc_offset(i)))*c_lunit
      enddo
     
      xv(4:6,1:np_local) = xv(4:6,1:np_local)/sqrt(a) * c_vunit

#endif

      g_time = a

      g_mass(1:6) = 0.0d0
      g_mass(2) = real(c_munit,8)
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
      print*,mpi_npart
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      g_npartTotal(1:6) = 0
      g_nhighword(1:6) = 0
      if(rank == 0 ) print*,"total N",mpi_nparttotal
      g_npartTotal(2) = iand(mpi_nparttotal,2**32-1)
      if(rank == 0) print*, g_npartTotal
      g_nhighword(2) = ishft(mpi_nparttotal,-32)
      if(rank == 0) print*,g_nhighword

      if(mod(rank,1) == 0) then
         open(unit=21,file=trim(output),form='unformatted',status='replace')

         write(21) g_npart, g_mass, g_time, g_redshift, g_flag_sfr, g_flag_feedback, g_npartTotal, &
              g_flag_cooling, g_num_files, g_boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam, &
              g_flag_stellarage, g_flag_metals, g_nhighword, g_filler
         write(21) xv(1:3,1:np_local)
         write(21) xv(4:6,1:np_local)
#ifdef LONGIDS
         write(21) PID(1:np_local)
#else
         write(21) int(PID(1:np_local),4)
#endif
         close(21)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      ! if(mod(rank,2) /= 0) then
      !    open(unit=21,file=trim(output),form='unformatted',status='replace')
      !    write(21) g_npart, g_mass, g_time, g_redshift, g_flag_sfr, g_flag_feedback, g_npartTotal, &
      !         g_flag_cooling, g_num_files, g_boxsize, g_Omega0, g_OmegaLambda, g_HubbleParam, &
      !         g_flag_stellarage, g_flag_metals, g_nhighword, g_filler
      !    write(21) xv(1:3,1:np_local)
      !    write(21) xv(4:6,1:np_local)
      !    write(21) PID(1:np_local)
      !    close(21)
      ! endif

      ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      deallocate(xv,PID)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine cube2gadget_call
end module cube2gadget_module

program main
  use mympi
  use cube2gadget_module
  implicit none
  integer :: i
  integer :: startsnap,stopsnap
  character(len=100) :: snapstr
  character(len=1000) :: prefix
  integer :: snap
  character(len=1000) :: input_dir,output_dir  
  call mympi_init
  call getarg(1,snapstr)
  read (snapstr, '(i10)') snap
  if(rank == 0) then
    i = 1
    open(22,file="./halofinds",action="read",status='old')
    open(23,file="./snap.txt",action="write",status="replace")
    redshift_list(1:maxsnap) = -1.0
    totalsnaps=0
100  read(22,fmt=*,end=200) redshift_list(i)
    write(23,*) 1./(1.+redshift_list(i))
    print*,redshift_list(i)
    totalsnaps = totalsnaps+1
    i = i+1
    goto 100
200  close(22)
    close(23)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call mpi_bcast(redshift_list,maxsnap,mpi_real,0,mpi_comm_world,ierr)
  call mpi_bcast(totalsnaps,1,mpi_integer,0,mpi_comm_world,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#include "config.f90"
  if(rank == 0) print*,"totalsnaps",totalsnaps
#ifndef IC
  call cube2gadget_call(input_dir,output_dir,prefix,snap)
#else
  call cube2gadget_call(input_dir,output_dir,prefix,-1)
#endif
  call mympi_end
end program main  







!
!  Copyright 2017 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!This file is "cmd_sc.f90"
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module control_sc_cmd
  use Global_Variables
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
  use salmon_parallel, only: nproc_id_global
  implicit none
contains

  subroutine init_cmd_sc
  use salmon_global
  use opt_variables
  use salmon_file
  use misc_routines
  use timer
  use force_field_CRK, only: allocate_ff_CRK,load_force_field_CRK,prepare_force_field_CRK
!  use inputoutput, only: au_length_aa
  implicit none
  integer :: is,ispecies, imol,num_mol
  integer :: tmp_natom, tmp_mol2species(NI), tmp_mol2atom_top(NI+1)
  character(1024) :: ifile_cmd_sc
  character(100)  :: ctmp1

  !(read input no.1)
  if (comm_is_root(nproc_id_global)) then
     write(*,*) "--- Read cmd_sc.inp.inp file for classical MD calculation"
     ifile_cmd_sc = "./cmd_sc.inp"
     open(800,file=trim(ifile_cmd_sc),status="old")
     read(800,*) !comment line
     read(800,*) force_field_system  ! force field of whole system
     read(800,*) !comment line
     read(800,*) nmol_s        ! # of molecular species in the system

     !check force field name
     if( force_field_system=='CRK') then
        write(*,*) "   force field is chosen to CRK"
     else 
        call Err_finalize("wrong force field name")
     endif

  endif
  call comm_bcast(force_field_system, nproc_group_global)
  call comm_bcast(nmol_s,             nproc_group_global)

  allocate( name_mol_s(nmol_s) ) ! molecular species name

  !(read input no.2)
  if (comm_is_root(nproc_id_global)) then
     !read molecular species and its name
     do is=1,nmol_s
        read(800,*) ispecies, ctmp1
        name_mol_s(ispecies) = ctmp1    ! molecular name e.x. "WAT"
     enddo
  endif
  call comm_bcast(name_mol_s, nproc_group_global)

  if(force_field_system=='CRK') then
     call allocate_ff_CRK
     do is = 1,nmol_s
        call load_force_field_CRK(is,name_mol_s(is))
     enddo
  endif

  !(read input no.3 & close the file)
  if (comm_is_root(nproc_id_global)) then
     read(800,*) !comment line
     tmp_natom = 0
     nmol  = 0
     tmp_mol2atom_top(1) = 1
     do 
        read(800,*,end=100) ispecies, num_mol
        tmp_mol2species(nmol+1:nmol+num_mol) = ispecies
        do imol=1,num_mol
           tmp_mol2atom_top(nmol+imol+1) = tmp_mol2atom_top(nmol+imol) + natom_mol_s(ispecies)
        enddo
        nmol = nmol + num_mol

        tmp_natom = tmp_natom + num_mol * natom_mol_s(ispecies)
        if(tmp_natom==NI) exit
        if(tmp_natom.gt.NI) then
           call Err_finalize("wrong in force field file")
           stop
        endif
     enddo
100  close(800)
  endif

  call comm_bcast(nmol, nproc_group_global)
  allocate( mol2species(nmol) )
  allocate( mol2atom_top(nmol))
  allocate( mol2atom_cnt(nmol))
  allocate( natom_mol(nmol) )
  if (comm_is_root(nproc_id_global)) mol2species(1:nmol) = tmp_mol2species(1:nmol)
  call comm_bcast(mol2species, nproc_group_global)
  if (comm_is_root(nproc_id_global)) mol2atom_top(1:nmol)= tmp_mol2atom_top(1:nmol)
  call comm_bcast(mol2atom_top, nproc_group_global)

  do imol=1,nmol
     natom_mol(imol) = natom_mol_s(mol2species(imol))
  enddo

  !(log)
  if (comm_is_root(nproc_id_global)) then
     write(*,*) "   total number of molecules =", nmol
     write(*,*) "   total number of molecular species =", nmol_s
     write(*,*) "   molecular species:"
     do is=1,nmol_s
        write(*,*) "  ",is, name_mol_s(is), natom_mol_s(is) 
     enddo
     write(*,*) "   molecules:"
     do imol=1,nmol
        write(*,*) "  ",imol,mol2species(imol), natom_mol(imol), mol2atom_top(imol)
     enddo
     write(*,*) "----------"
  endif

  call prepare_force_field_CRK

  end subroutine

!------------------------------------------------------------------------
  subroutine opt_cmd
    use md_ground_state
    use optimization, only: cal_mean_max_forces
    use force_field_CRK, only: cal_force_energy_CRK
    implicit none
    integer :: iter, ia
    real(8) :: fmax_conv, fmax, fave, step_size, opt_direction(3,NI), zero
    character(100) :: comment_line

    fmax_conv = convrg_opt_fmax
!    step_size = 0.002d0  !test
    step_size = 0.0005d0  !test
    zero = 0d0

    call cal_force_energy_CRK
    call cal_mean_max_forces(NI,force,fave,fmax)
    opt_direction(:,:) = force(:,:)/fave

    ! Export to file_trj (initial step)
    if (out_rvf_rt=='y')then
       write(comment_line,110) -1
       call write_xyz(comment_line,"new","rvf")
       call write_xyz(comment_line,"add","rvf")
    endif


    if(comm_is_root(nproc_id_global)) then
       write(*,'(a)') "# step, Uene, fmax, fave [in a.u. ]"
    endif

    call comm_sync_all


    ! === Main loop of optimization step ===
    do iter = 0, Nt

       do ia=1,NI
          Rion(:,ia) = Rion(:,ia) + opt_direction(:,ia)*step_size
       enddo

       call cal_force_energy_CRK
       call cal_mean_max_forces(NI,force,fave,fmax)
       opt_direction(:,:) = force(:,:)/fave

       !check convergence
       if(abs(fmax).lt.fmax_conv) then
          if(comm_is_root(nproc_id_global)) then
             write(*,*) "Optimization in CMD was Converged: fmax=", fmax
          endif
          exit
       endif

       !---Write section---
       ! Export to file_rt_data
       call write_cmd_data(iter,Uene,fmax,fave,zero,zero,zero,"opt")

       ! Export to standard log file
       if(comm_is_root(nproc_id_global)) then
          write(*,120) iter ,Uene, fmax, fave
120       format(1x,i6, 100e20.10E3)
       endif

       ! Export to file_trj
       if (out_rvf_rt=='y' .and. mod(iter,out_rvf_rt_step)==0)then
          write(comment_line,110) iter
110       format("#cmd  step=",i8)
          if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
               &  write(comment_line,112) trim(comment_line), xi_nh
112       format(a,"  xi_nh=",e18.10)
          call write_xyz(comment_line,"add","rvf")
       endif

    enddo  !main loop of it

    call comm_sync_all

    if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'


  end subroutine

!------------------------------------------------------------------------
  subroutine cmd_sc
    use md_ground_state
    use force_field_CRK, only: cal_force_energy_CRK
    implicit none
    integer :: it,ia
    real(8) :: dt_h, aforce(3,NI)
    real(8) :: Temp_sys, mass_au
    real(8) :: Htot, Enh, Enh_gkTlns, gkT, Qnh
    character(100) :: comment_line

    call take_back_mol_into_ucell

    ! Export to file_trj (initial step)
    if (out_rvf_rt=='y')then
       call cal_force_energy_CRK
       write(comment_line,110) -1, 0.0d0
       if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
       &  write(comment_line,112) trim(comment_line), xi_nh
       call write_xyz(comment_line,"new","rvf")
       call write_xyz(comment_line,"add","rvf")
    endif

    it         = 0
    dt_h       = dt*0.5d0
    Enh_gkTlns = 0d0
    Enh        = 0d0

    if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
       gkT = 3d0*NI * kB/hartree2J*temperature0_ion
       Qnh = gkT * thermostat_tau**2d0
    endif

    if(comm_is_root(nproc_id_global)) then
       write(*,'(a)') "# time, Tene, Uene, Eall, Enh, Htot, Temperature [in a.u. & K]"
    endif

    call comm_sync_all


    ! === Main loop of time step ===
    do it = 0, Nt

       ! Velocity Verlet integrator

       !NHC act on velocity with dt/2
       if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
          call apply_nose_hoover_velocity(dt_h)
       endif

       do ia=1,NI
          mass_au = umass*Mass(Kion(ia))
          velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       enddo

       if(step_velocity_scaling>=1 .and. mod(it,step_velocity_scaling)==0) then
          call cal_Tion_Temperature_ion(Tene,Temp_sys,velocity)
          call apply_velocity_scaling_ion(Temp_sys,velocity)
       endif

       do ia=1,NI
          Rion(:,ia) = Rion(:,ia) + velocity(:,ia)*dt
       enddo
       call take_back_mol_into_ucell

       !put SHAKE here in future if needed

       !NHC act on thermostat with dt
       if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
          call cal_Tion_Temperature_ion(Tene,Temp_sys,velocity)
          call apply_nose_hoover_thermostat(Temp_sys,dt)
          Enh_gkTlns = Enh_gkTlns + gkT * xi_nh*dt
          Enh        = Enh_gkTlns + 0.5d0 * Qnh * xi_nh*xi_nh
       endif


       !update force (electric state) with updated coordinate
       aforce(:,:) = force(:,:)
       call cal_force_energy_CRK

       aforce(:,:) = 0.5d0*( aforce(:,:) + force(:,:) )

       do ia=1,NI
          mass_au = umass*Mass(Kion(ia))
          velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       enddo

       !NHC act on velocity with dt/2
       if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
          call apply_nose_hoover_velocity(dt_h)
       endif

       !remove system momentum
       if(stop_system_momt=='y') call remove_system_momentum(0)


       call cal_Tion_Temperature_ion(Tene,Temp_sys,velocity)
       Eall     = Uene + Tene
       Htot     = Eall + Enh


       !---Write section---
       ! Export to file_rt_data
       call write_cmd_data(it,Tene,Uene,Eall,Temp_sys,Enh,Htot,"rt")

       ! Export to standard log file
       if(comm_is_root(nproc_id_global)) then
          write(*,120) it*dt,Tene,Uene,Eall,Enh,Htot,Temp_sys
120       format(1x,f10.4, 7e20.10E3,f12.3)
       endif

       ! Export to file_trj
       if (out_rvf_rt=='y' .and. mod(it,out_rvf_rt_step)==0)then
          write(comment_line,110) it, it*dt
110       format("#cmd  step=",i8,"   time",e16.6)
          if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
               &  write(comment_line,112) trim(comment_line), xi_nh
112       format(a,"  xi_nh=",e18.10)
          call write_xyz(comment_line,"add","rvf")
       endif


    enddo  !main loop of it

    call print_restart_data_md_gs
    call comm_sync_all

    if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'

  end subroutine

  subroutine take_back_mol_into_ucell
    implicit none
    integer j,imol,iatom_c,iatom,ia

    do imol=1,nmol
       iatom_c = mol2atom_cnt(imol)
       do j=1,3
          if( Rion(j,iatom_c).lt.0d0 ) then
             do ia= 1,natom_mol(imol)
                iatom = mol2atom_top(imol) + ia-1
                Rion(j,iatom) = Rion(j,iatom) + aL(j)
             enddo
          else if( Rion(j,iatom_c).gt.aL(j) ) then
             do ia= 1,natom_mol(imol)
                iatom = mol2atom_top(imol) + ia-1
                Rion(j,iatom) = Rion(j,iatom) - aL(j)
             enddo
          endif
       enddo
    enddo

  end subroutine

  subroutine write_cmd_data(it,dat1,dat2,dat3,dat4,dat5,dat6,mode)
    use inputoutput, only: t_unit_time,t_unit_energy
    implicit none
    integer :: fh_rt=202, it
    real(8) :: dat1,dat2,dat3,dat4,dat5,dat6
    real(8) :: Tene_t,Uene_t,Eall_t,Temp_sys_t,Enh_t,Htot_t, Fmax_t,Fave_t
    character(*) :: mode

    if(comm_is_root(nproc_id_global)) then

       if(mode=="rt")then         
          Tene_t     = dat1
          Uene_t     = dat2
          Eall_t     = dat3
          Temp_sys_t = dat4
          Enh_t      = dat5
          Htot_t     = dat6
       else if(mode=="opt")then
          Uene_t = dat1
          Fmax_t = dat2
          Fave_t = dat3
       endif

      if(it==0) then
         open(fh_rt,file=trim(file_rt_data),status="unknown")

         if(mode=="rt")then         
            write(fh_rt, '("#",1X,A)') "time (Classical MD)"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Tene", "Kinetic energy"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Uene", "Potential energy"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Eall", "Total energy"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Temp_sys", "Temperature of system"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Enh", "Energy of NH thermostat"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Hnvt", "Hamiltonian with NH thermostat"
         
            write(fh_rt, '("#",99(1X,I0,":",A,"[",A,"]"))')&
                 &  1, "Time",   trim(t_unit_time%name),   &
                 &  2, "Tene",   trim(t_unit_energy%name), &
                 &  3, "Uene",   trim(t_unit_energy%name), &
                 &  4, "Eall",   trim(t_unit_energy%name), &
                 &  5, "Temp_sys", "K",                    &
                 &  6, "Enh",    trim(t_unit_energy%name), &
                 &  7, "Hnvt",   trim(t_unit_energy%name)

         else if(mode=="opt")then
            write(fh_rt, '("#",1X,A)') "step (Optimization)"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Uene", "Potential energy"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Fmax", "Maxmum Force on atom"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Fave", "Average Force on atom"

            write(fh_rt, '("#",99(1X,I0,":",A,"[",A,"]"))')&
                 &  1, "Step",   "",                       &
                 &  2, "Uene",   trim(t_unit_energy%name), &
                 &  3, "Fmax",   "a.u.",                   &
                 &  4, "Fave",   "a.u."
         endif

      endif

      if(mode=="rt")then         
         write(fh_rt, "(F16.8,99(1X,E23.15E3))") &
              & it * dt   * t_unit_time%conv,   &
              & Tene_t    * t_unit_energy%conv, &
              & Uene_t    * t_unit_energy%conv, &
              & Eall_t    * t_unit_energy%conv, &
              & Temp_sys_t,                     &
              & Enh_t     * t_unit_energy%conv, &
              & Htot_t    * t_unit_energy%conv
      else if(mode=="opt")then
         write(fh_rt, "(I8,99(1X,E23.15E3))")   &
              & it,                             &
              & Uene_t    * t_unit_energy%conv, &
              & Fmax_t,                         &
              & Fave_t
      endif

      flush(fh_rt)

      if(it==Nt) close(fh_rt)
    end if

    call comm_sync_all
    return
  end subroutine write_cmd_data


end module

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
  use force_field_CRK, only: allocate_ff_CRK,load_force_field_CRK,prepare_force_field_CRK,&
                              iter_max, dQave_thresh, a_ewald, r_cutoff, G_cutoff
  use inputoutput, only: au_length_aa
  implicit none
 !logical :: flag_advanced
  integer :: is,ispecies, imol,num_mol
  integer :: tmp_natom, tmp_mol2species(NI), tmp_mol2atom_top(NI+1)
  character(1024) :: ifile_cmd_sc
  character(100)  :: ctmp1 !keywd

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

  !(read input no.3)
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
100  continue
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

  !(read input no.4)
  if (comm_is_root(nproc_id_global)) then
     read(800,*) iter_max       !=50 is enough
     read(800,*) dQave_thresh   !=1d-7 [au]  (from paper)
     read(800,*) a_ewald        != [1/A]
     read(800,*) r_cutoff       !=13.0 [A]   (from paper)
     read(800,*) G_cutoff       !=1.47 [1/A] (from paper)
     a_ewald  = a_ewald *au_length_aa !to [1/Bohr]
     r_cutoff = r_cutoff/au_length_aa !to [Bohr]
     G_cutoff = G_cutoff*au_length_aa !to [1/Bohr]
  endif

  call comm_bcast(iter_max,     nproc_group_global)
  call comm_bcast(dQave_thresh, nproc_group_global)
  call comm_bcast(a_ewald,      nproc_group_global)
  call comm_bcast(r_cutoff,     nproc_group_global)
  call comm_bcast(G_cutoff,     nproc_group_global)

  !(read input no.5 & close the file)
  if (comm_is_root(nproc_id_global)) then
!     read(800,*) flag_advanced  !=.false.
!     if(flag_advanced)then
!        read(800,*) keywd
!     endif
     close(800)
  endif

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
    implicit none

    call  opt_cmd_conjugate
!    call  opt_cmd_steepest

  end subroutine
!------------------------------------------------------------------------
  subroutine opt_cmd_conjugate   !conjugated method
    use md_ground_state
    use optimization, only: variable_3xNto3N,cal_inner_product, &
                      get_predicted_zero,update_2pt_line_search,cal_mean_max_forces
    use force_field_CRK, only: cal_force_energy_CRK
    implicit none
    logical :: flag_accept
    integer :: i,ia, iter,iter_perp,iter_line, Nopt_perp,Nopt_line
    real(8) :: gm, tmp1,tmp2,tmp3, Rion_save(3,NI),dRion_rmsd
    real(8) :: StepLen_line(2), StepLen_line_small
    real(8) :: StepLen_line_new, StepLen_line_zero, sl_zero
    real(8) :: StepLen_line_Up, StepLen_line_Dw
    real(8) :: SearchDirection(3,NI),SearchDirection_1d(3*NI)
    real(8) :: Uene_prev,Uene_prev_line,Uene_zero,Uene_new,dUene
    real(8) :: fmax_conv,fmax,fave, zero,zero3(3)
    real(8) :: force_prev(3,NI), force_1d(3*NI), force_prev_1d(3*NI)
    real(8) :: F_line_conv, F_line_2pt(2), F_line,F_line_new,F_line_zero
    character(100) :: comment_line

    zero     = 0d0
    zero3(:) = 0d0

!    Nopt_perp = 1000
    Nopt_perp = Nt
    Nopt_line = 100
    StepLen_line_small = 0.1d0   !(small step to guess initial step length)
    StepLen_line_Up    = cg_alpha_up
    StepLen_line_Dw    = cg_alpha_down
    F_line_conv        = convrg_opt_fmax
    fmax_conv          = convrg_opt_fmax
    
    if(comm_is_root(nproc_id_global)) then
       write(*,*) "  [Set following in optimization]"
       write(*,*) "  Max optimization CG step    =",Nopt_perp
       write(*,*) "  Max line search opt step    =",Nopt_line
       write(*,*) "  Up rate of line-search step length  =",real(StepLen_line_Up)
       write(*,*) "  Down rate of line-search step length=",real(StepLen_line_Dw)
       write(*,*) "  Convergence threshold of F_line =",real(F_line_conv)
       write(*,*) "  Convergence threshold of Fmax   =",real(fmax_conv)
    endif
    
    !Initial Step Procedure
    call cal_force_energy_CRK
    SearchDirection(:,:) = force(:,:)
    do ia=1,NI
       if(flag_geo_opt_atom(ia)=='n') SearchDirection(:,ia)=0d0  !fix atom
    enddo
    call variable_3xNto3N(NI,force,force_1d)
    call variable_3xNto3N(NI,SearchDirection,SearchDirection_1d)
    call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)

    write(comment_line,110) 0, 0
    call write_xyz(comment_line,"new","r  ")
    call cal_mean_max_forces(NI,force,fave,fmax)
    if(comm_is_root(nproc_id_global)) write(*,135) 0,0,fmax,fave,Uene,0d0

    !(initial check of convergence: if force is enough small, exit before opt)
    if(fmax .le. fmax_conv) then
       if(comm_is_root(nproc_id_global)) &
            &  write(*,*) " Max force is enough small: stop calculation"
       call Err_finalize(" Max force is enough small: stop calculation")
       stop
    endif

    !==== Main Loop ====
    !---iteration for perpendicular direction---
    do iter_perp =1,Nopt_perp   

       if(comm_is_root(nproc_id_global)) then
          write(*,*) "==================================================="
          write(*,*) "CG Optimization Step = ", iter_perp
       endif

       !previous value
       force_prev(:,:) = force(:,:)
       Uene_prev = Uene
       Uene_prev_line = Uene

       !(store)
       Rion_save(:,:)= Rion(:,:)

       !Set initial region to be searched (=> set initial two points)
       !(calculate forces at 2 points and adjust initial step length)
       call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
       StepLen_line(1) = 0d0
       F_line_2pt(1)   = F_line

       if(F_line_2pt(1).le.0d0) then
          if(comm_is_root(nproc_id_global)) &
          &  write(*,*) "WARNING: Search Direction is opposite"
          exit
       endif

       !(guess initial step length from numerical gradient of force)
       i= 2
       StepLen_line(i) = StepLen_line_small
       do
          Rion(:,:)   = Rion_save(:,:) + StepLen_line(i)* SearchDirection(:,:)
          Rion_eq(:,:)= Rion(:,:)
          call cal_force_energy_CRK
          call variable_3xNto3N(NI,force,force_1d)
          call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
          F_line_2pt(i) = F_line
          call get_predicted_zero(StepLen_line,F_line_2pt,sl_zero)
          dRion_rmsd  = sl_zero*sqrt(sum(SearchDirection_1d(:)**2)/NI)
          if(comm_is_root(nproc_id_global)) &
          & write(*,104)" guessed-initial-alpha=",sl_zero,"  dRion-RMSD=",dRion_rmsd

          if(abs(dRion_rmsd).lt.1d-6 .or. dRion_rmsd.gt.1d0) then
             ! just ad-hoc treatment (unusual case: I don't know why this case exists)
             StepLen_line(i) = StepLen_line_small*5d0
             if(comm_is_root(nproc_id_global)) then
                write(*,*) " Could not find good guess of initial alpha: too small or large RMSD"
                write(*,*) " --> put guessed initial alpha=",real(StepLen_line(i))
             endif
             exit
          endif

          if(sl_zero.lt.0d0) then
             StepLen_line(i) = StepLen_line(i) * StepLen_line_Dw
          else 
             StepLen_line(i) = sl_zero
             exit
          endif
       enddo

       !(adjust initial step length)
10     continue
       i= 2
       Rion(:,:)   = Rion_save(:,:) + StepLen_line(i)* SearchDirection(:,:)
       Rion_eq(:,:)= Rion(:,:)
       dRion_rmsd  = StepLen_line(i)*sqrt(sum(SearchDirection_1d(:)**2)/NI)
       call cal_force_energy_CRK
       call variable_3xNto3N(NI,force,force_1d)
       call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
       F_line_2pt(i) = F_line
       
       !(adjust initial step length so as to include zero-force between 2 points)
       if( F_line_2pt(1)*F_line_2pt(2).gt.0d0 ) then
          !Here, this algorithm can be replaced by linearly predicted point(step length)
          !if a lot of "up", modify here with prediction
          StepLen_line(1) = StepLen_line(2)
          F_line_2pt(1)   = F_line_2pt(2)
          StepLen_line(2) = StepLen_line(2) * StepLen_line_Up
          dRion_rmsd      = StepLen_line(2)*sqrt(sum(SearchDirection_1d(:)**2)/NI)
          if(comm_is_root(nproc_id_global)) &
          & write(*,104)" adjusting initial alpha(up)->",StepLen_line(2),"  dRion-RMSD=",dRion_rmsd
          goto 10
       else
          if(comm_is_root(nproc_id_global)) &
          & write(*,104)" initial alpha=",StepLen_line(2),"  dRion-RMSD=",dRion_rmsd
       endif


       !--- iteration loop for line search to find zero-force ---
       do iter_line= 1,Nopt_line
          
          !(narrow search range)
          call get_predicted_zero(StepLen_line,F_line_2pt,sl_zero)

          if(abs(F_line_2pt(1)/F_line_2pt(2)).gt.1d0) then
             StepLen_line_new= 0.5d0*( StepLen_line(1) + sl_zero )
          else
             StepLen_line_new= StepLen_line(1) + 0.5d0*(StepLen_line(2)-sl_zero)
          endif

          !(calculate electronic state at the new point)
          Rion(:,:)= Rion_save(:,:) + StepLen_line_new* SearchDirection(:,:)
          Rion_eq(:,:)= Rion(:,:)
          write(comment_line,110) iter_perp, iter_line
          call write_xyz(comment_line,"add","r  ")
          call cal_force_energy_CRK
          call variable_3xNto3N(NI,force,force_1d)
          call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
          F_line_new = F_line
          Uene_new   = Uene

          flag_accept=.false.
          if(abs(F_line_2pt(1)/F_line_2pt(2)).gt.1d0) then
             if( F_line_new*F_line_2pt(1) .gt. 0d0 ) flag_accept=.true.
          else
             if( F_line_new*F_line_2pt(2) .gt. 0d0 ) flag_accept=.true.
          endif
          
          if(flag_accept) &
          &  call update_2pt_line_search(StepLen_line,F_line_2pt,StepLen_line_new,F_line_new)


          !(get predicted zero-crossin coordinate)
          call get_predicted_zero(StepLen_line,F_line_2pt,StepLen_line_zero)

          !(calculate electronic state at the predicted zero)
          Rion(:,:)   = Rion_save(:,:) + StepLen_line_zero* SearchDirection(:,:)
          Rion_eq(:,:)= Rion(:,:)
          dRion_rmsd  = StepLen_line_zero*sqrt(sum(SearchDirection_1d(:)**2)/NI)
          write(comment_line,110) iter_perp, iter_line
          call write_xyz(comment_line,"add","r  ")
          call cal_force_energy_CRK
          call variable_3xNto3N(NI,force,force_1d)
          call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
          F_line_zero = F_line
          Uene_zero   = Uene
          
          !(log)
          if(comm_is_root(nproc_id_global)) then
             write(*,100) "   alpha-line:2pt=",(StepLen_line(i),i=1,2),&
                        & " |predict-zero=",    StepLen_line_zero
             write(*,100) "   force-line:2pt=",(F_line_2pt(i),i=1,2), &
                        & " |predict-zero=",    F_line_zero
          endif

          !Judge Convergence for line search opt
          dUene   = Uene_zero - Uene_prev_line
          if(comm_is_root(nproc_id_global)) &
          & write(*,130)iter_perp,iter_line,StepLen_line_zero,F_line_zero,dUene,dRion_RMSD

          if(abs(F_line_zero) .le. F_line_conv)then
             if(comm_is_root(nproc_id_global)) then
                write(*,*)
                write(*,*) "Converged(line search opt) in perpendicular-step",iter_perp
             endif
             exit
          else if(iter_line==Nopt_line) then
             if(comm_is_root(nproc_id_global)) then
                write(*,*) "Not Converged(line search opt) in perpendicular-step",iter_perp
                write(*,*) "==================================================="
             endif
             exit
          endif

          !(update two points among current predicted-zero and two points)
          call update_2pt_line_search(StepLen_line,F_line_2pt,StepLen_line_zero,F_line_zero)

          !(preparation for next cycle)
          Uene_prev_line = Uene_zero

       enddo

       !Force calculation
       call cal_force_energy_CRK
       call variable_3xNto3N(NI,force,force_1d)
       call cal_mean_max_forces(NI,force,fave,fmax)

       !Judge Convergence for perpendicular opt (main judge)
       dUene = Uene - Uene_prev
       if(comm_is_root(nproc_id_global)) &
       &  write(*,135) iter_perp, iter_line, fmax,fave,Uene,dUene
       if(abs(Fmax) .le. fmax_conv) then
          if(comm_is_root(nproc_id_global)) then
             write(*,*) "==================================================="
             write(*,*)
             write(*,*) "Optimization Converged"
             write(*,150) iter_perp,fmax,fave,Uene
             write(*,*) "==================================================="
          endif
          exit
       else if(iter_perp==Nopt_perp) then
          if(comm_is_root(nproc_id_global)) then
             write(*,*) "==================================================="
             write(*,*)
             write(*,*) "Optimization Did Not Converged"
             write(*,150) iter_perp,fmax,fave,Uene
             write(*,*) "==================================================="
          endif
       endif
       

       !Update search direction vector for perpendicular step
       call variable_3xNto3N(NI,force,force_1d)
       call variable_3xNto3N(NI,force_prev,force_prev_1d)
       call cal_inner_product(3*NI,force_1d,force_1d,tmp1)
       call cal_inner_product(3*NI,force_prev_1d,force_prev_1d,tmp2)
       call cal_inner_product(3*NI,force_1d,force_prev_1d,tmp3)
       !gm = tmp1/tmp2         !(by Fletcher-Reeves)
       gm = (tmp1-tmp3)/tmp2   !(by Polak-Ribiere)--usually best, but sometimes direction is not downward.
       SearchDirection(:,:) = force(:,:) + gm * SearchDirection(:,:)
       do ia=1,NI
          if(flag_geo_opt_atom(ia)=='n') SearchDirection(:,ia)=0d0  !fix atom
       enddo
       call variable_3xNto3N(NI,SearchDirection,SearchDirection_1d)

       !---Write section---
       ! Export to file_rt_data
       call write_cmd_data(iter,Uene,fmax,fave,zero,zero,zero,zero,zero3,zero3,zero3,"opt")


    enddo !end of opt iteraction========================

    velocity(:,:) = 0d0
    call print_restart_data_md_gs
    call comm_sync_all

    if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'


100 format(a,2e13.5,a,e13.5)
104 format(a,e13.5,a,e13.5)
110 format("#opt   step-perp=",i4,"   step-line",i4)
120 format(a,e14.6,a,e14.6)
130 format(" step=",i3," -",i3,"  alpha=",e12.4,"  F-line=",e12.4,"  dE-line=",e12.4,"  dRion-RMSD=",e12.4)
135 format(" step=(perp=",i3,", line=",i3,")  Fmax=",e11.4,"  Fave=",e11.4,"  E=",e16.8,"  dE-perp=",e12.4)
150 format("    Iteration step=",i5,    / &
    &      "    Force(maximum)=",e18.10," [a.u.]",/ &
    &      "    Force(average)=",e18.10," [a.u.]",/ &
    &      "    Total Energy  =",e18.10," [a.u.]" )

  end subroutine

!------------------------------------------------------------------------
  subroutine opt_cmd_steepest  !not used now
    use md_ground_state
    use optimization, only: cal_mean_max_forces
    use force_field_CRK, only: cal_force_energy_CRK
    implicit none
    integer :: iter, ia
    real(8) :: fmax_conv,fmax,fave, step_size, opt_direction(3,NI),zero,zero3(3)
    character(100) :: comment_line

    fmax_conv = convrg_opt_fmax
!    step_size = 0.002d0  !test
!    step_size = 0.001d0  !test
    step_size = 0.0005d0  !test
    zero     = 0d0
    zero3(:) = 0d0

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
       call write_cmd_data(iter,Uene,fmax,fave,zero,zero,zero,zero,zero3,zero3,zero3,"opt")

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

    velocity(:,:) = 0d0
    call print_restart_data_md_gs
    call comm_sync_all

    if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'

  end subroutine

!------------------------------------------------------------------------
  subroutine cmd_sc
    use md_ground_state
    use force_field_CRK, only: cal_force_energy_CRK, Vuc
    implicit none
    integer :: it,ia
    real(8) :: dt_h, aforce(3,NI)
    real(8) :: Temp_sys, mass_au
    real(8) :: Htot, Enh, Enh_gkTlns, gkT, Qnh, Eem
    real(8) :: Mt(3),Jmt(3)
    character(100) :: comment_line

    call take_back_mol_into_ucell
    call init_Ac
    it = 0
    Et(:) = -( Ac_ext(it+1,:) - Ac_ext(it-1,:) ) / (2*dt)

    ! Export to file_trj (initial step)
    if (out_rvf_rt=='y')then
       call cal_force_energy_CRK
       write(comment_line,110) -1, 0.0d0
       if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
       &  write(comment_line,112) trim(comment_line), xi_nh
       call write_xyz(comment_line,"new","rvf")
       call write_xyz(comment_line,"add","rvf")
    endif

    Mt(:)=0d0
    call cal_dipole_moment_current(Mt,Jmt,dt)
    Jmt(:)=0d0

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

       !Electric field
       Et(:) = -( Ac_ext(it+1,:) - Ac_ext(it-1,:) ) / (2*dt)
       Eem   = Vuc * sum(Et(:)**2)/(8d0*pi)

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
       Eall     = Uene + Tene + Eem
       Htot     = Eall + Enh


       !---Analyses section---
       call cal_dipole_moment_current(Mt,Jmt,dt)

       !---Write section---
       ! Export to file_rt_data
       call write_cmd_data(it,Tene,Uene,Eem,Eall,Temp_sys,Enh,Htot,Et,Mt,Jmt,"rt")

       ! Export to standard log file
       if(comm_is_root(nproc_id_global)) then
          write(*,120) it*dt,Tene,Uene,Eall,Enh,Htot,Temp_sys
120       format(1x,f14.4, 7e20.10E3,f12.3)
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

  subroutine write_cmd_data(it,dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,mode)
    use inputoutput, only: t_unit_time,t_unit_energy,t_unit_elec,t_unit_current
    implicit none
    integer :: fh_rt=202, it
    real(8) :: dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8(3),dat9(3),dat10(3)
    real(8) :: Tene_t,Uene_t,Eem_t,Eall_t,Temp_sys_t,Enh_t,Htot_t, Fmax_t,Fave_t
    real(8) :: E_t(3), M_t(3), Jm_t(3)
    character(*) :: mode

    if(comm_is_root(nproc_id_global)) then

       if(mode=="rt")then         
          Tene_t     = dat1
          Uene_t     = dat2
          Eem_t      = dat3
          Eall_t     = dat4
          Temp_sys_t = dat5
          Enh_t      = dat6
          Htot_t     = dat7
          E_t(:)     = dat8(:)
          M_t(:)     = dat9(:)
          Jm_t(:)    = dat10(:)
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
            write(fh_rt, '("#",1X,A,":",1X,A)') "Eem",  "Energy of electric field"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Eall", "Total energy"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Temp", "Temperature"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Enh",  "Energy of NH thermostat"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Hnvt", "Hamiltonian with NH thermostat"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Et", "Electric field"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Mt", "Dipole Moment of the system"
            write(fh_rt, '("#",1X,A,":",1X,A)') "Jmt","Matter current density"
         
            write(fh_rt, '("#",99(1X,I0,":",A,"[",A,"]"))')&
                 &  1, "Time",   trim(t_unit_time%name),   &
                 &  2, "Tene",   trim(t_unit_energy%name), &
                 &  3, "Uene",   trim(t_unit_energy%name), &
                 &  4, "Eem",    trim(t_unit_energy%name), &
                 &  5, "Eall",   trim(t_unit_energy%name), &
                 &  6, "Temp",   "K",                      &
                 &  7, "Enh",    trim(t_unit_energy%name), &
                 &  8, "Hnvt",   trim(t_unit_energy%name), &
                 &  9, "Et_x",   trim(t_unit_elec%name),   &
                 & 10, "Et_y",   trim(t_unit_elec%name),   &
                 & 11, "Et_z",   trim(t_unit_elec%name),   &
                 & 12, "Mt_x",   "a.u.",  &
                 & 13, "Mt_y",   "a.u.",  &
                 & 14, "Mt_z",   "a.u.",  &
                 & 15, "Jmt_x", trim(t_unit_current%name), &
                 & 16, "Jmt_y", trim(t_unit_current%name), &
                 & 17, "Jmt_z", trim(t_unit_current%name)

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
              & Eem_t     * t_unit_energy%conv, &
              & Eall_t    * t_unit_energy%conv, &
              & Temp_sys_t,                     &
              & Enh_t     * t_unit_energy%conv, &
              & Htot_t    * t_unit_energy%conv, &
              & E_t(1:3)  * t_unit_elec%conv,   &
              & M_t(1:3),                       &
              & Jm_t(1:3) * t_unit_current%conv
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

  subroutine cal_dipole_moment_current(dipole_mom,current_matter,dt)
    ! dipole_mom : (in & out)
    use force_field_CRK, only: Qai,Rxyz_wX,natom_wX_mol,Vuc
    implicit none
    integer :: imol,ia_wX
    real(8) :: dipole_mom(3), dipole_mom_last(3), current_matter(3), dt

    dipole_mom_last(:) = dipole_mom(:)

    !dipole moment in unit cell
    dipole_mom(:) = 0d0
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       dipole_mom(:) = dipole_mom(:) + Qai(ia_wX,imol)*Rxyz_wX(:,ia_wX,imol)
    enddo
    enddo

    !matter current defined by positive charge-->minus sign
    current_matter(:) = -(dipole_mom(:) - dipole_mom_last(:))/dt/Vuc


  end subroutine cal_dipole_moment_current

end module

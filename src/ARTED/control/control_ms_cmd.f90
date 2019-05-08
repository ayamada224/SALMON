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
!This file is "cmd_ms.f90"
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module control_ms_cmd
  use Global_Variables
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
  use salmon_parallel, only: nproc_id_global
  implicit none
  logical :: flag_ms_ff_LessPrint,  flag_fix_atoms_md 
  logical :: flag_use_absorb_bound, flag_use_light_source
  logical :: flag_method_newton_Et
  logical,allocatable :: flag_conv_iter2(:)
  integer :: iter_save, iter2, ix_light_source
  integer :: M_pml,nx_pml, nx1_pml,nx2_pml
  real(8),allocatable :: Rion_eq0(:,:), Rxyz_wX_ms(:,:,:,:)
  real(8),allocatable :: Qai_ms(:,:,:), Qai_old_ms(:,:,:)
  real(8),allocatable :: Qai_bk_ms(:,:,:)
  real(8),allocatable :: Uene_ms(:), Tene_ms(:), Eall_ms(:)
  real(8),allocatable :: Uene_intra_ms(:), force_intra_m(:,:,:)
  real(8),allocatable :: E_ms(:,:,:,:), B_ms(:,:,:,:)
  real(8),allocatable :: E_m(:,:), B_m(:,:), E_m_iter(:,:), E_m_bk(:,:)
  real(8),allocatable :: E_last(:,:,:,:), B_last(:,:,:,:)
  real(8),allocatable :: E_hist_m(:,:,:)
  real(8),allocatable :: D1_pml(:,:,:,:), D2_pml(:,:,:,:)
  real(8),allocatable :: D1_last(:,:,:,:),D2_last(:,:,:,:)
  real(8),allocatable :: velocity_m_mid(:,:,:),Rion_m_mid(:,:,:)
  real(8),allocatable :: Rxyz_wX_ms_mid(:,:,:,:), Rxyz_wX_ms_bk(:,:,:,:)
  real(8),allocatable :: data_local_MD(:,:,:)
  real(8) :: R_pml, sgm_max_E_pml, sgm_max_B_pml, a_pml
  character(3) :: ABC_method
  real(8),allocatable :: xi_nh_m(:),Enh_m_t(:,:),Hnvt_m_t(:,:)
  real(8) :: fric, conv_iter2, dEt_pc
  integer :: npt_intpl_E
  real(8),allocatable :: vpr(:,:,:), vmr(:,:,:), rvec(:,:,:), Qvec_bk(:,:)
  real(8),allocatable :: vec2_tmp_m(:,:),vec3_tmp_m(:,:)
contains

subroutine init_cmd_maxwell_ms
  use salmon_global
  use Global_Variables
  use opt_variables
  use salmon_parallel
  use salmon_communication
  use salmon_file
  use misc_routines
  use timer
  use force_field_CRK, only: allocate_ff_CRK,load_force_field_CRK,prepare_force_field_CRK,&
                              iter_max, a_ewald, r_cutoff, G_cutoff, NI_wX  !, dQave_thresh
  use inputoutput, only: au_length_aa, au_time_fs
  implicit none
  integer :: is,ispecies, imol,num_mol, ia,j,ik
  integer :: tmp_natom, tmp_mol2species(NI), tmp_mol2atom_top(NI+1)
  real(8) :: tmpr(3), uconv, tau_fric
  character(1024) :: ifile_cmd_ms
  character(100)  :: char_atom, ctmp1
  
  !keyword check
  if(use_ms_maxwell=='n') then
     call end_parallel
     stop
  endif

  flag_method_newton_Et = .true.  !must be .true.

  npt_intpl_E = 4     !number of interpolation for E(t) to predictor
 !(read from input now)
 !conv_iter2  = 1d-10 !convergence of E(t) for pred-corrector [au]
  if(comm_is_root(nproc_id_global)) then
     write(*,*) "number of interpolation points for E(t) prediction=", npt_intpl_E
     write(*,*) "convergence for predictor-corrector of E(t) [au]=", conv_iter2
  endif
  allocate( E_ms(1:3,mx1_m:mx2_m,my1_m:my2_m,mz1_m:mz2_m) )
  allocate( B_ms(1:3,mx1_m:mx2_m,my1_m:my2_m,mz1_m:mz2_m) )
  allocate( E_m(1:3,1:nmacro) )
  allocate( B_m(1:3,1:nmacro) )
  allocate( E_last(1:3,mx1_m:mx2_m,my1_m:my2_m,mz1_m:mz2_m) )
  allocate( B_last(1:3,mx1_m:mx2_m,my1_m:my2_m,mz1_m:mz2_m) )
  allocate( E_hist_m(1:3,1:nmacro,npt_intpl_E-1) )
  allocate( E_m_bk(1:3,1:nmacro), E_m_iter(1:3,1:nmacro) )
  allocate( flag_conv_iter2(1:nmacro) )
  E_hist_m(:,:,:)=0d0
  allocate( velocity_m_mid(3,NI,nmacro_s:nmacro_e) )
  allocate( Rion_m_mid(3,NI,nmacro_s:nmacro_e) )

  allocate( data_local_MD(3,nmacro_s:nmacro_e,0:Nt) )

  if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
     allocate( Enh_m_t(0:Nt,nmacro_s:nmacro_e) )
     allocate( Hnvt_m_t(0:Nt,nmacro_s:nmacro_e) )
     allocate( xi_nh_m(nmacro_s:nmacro_e) )
     xi_nh_m(:) = 0d0   !always start from zero even if restart
  endif

  !(read input no.1)
  if (comm_is_root(nproc_id_global)) then
     write(*,*) "--- Read cmd_ms.inp.inp file for classical MD multiscale calculation"
     ifile_cmd_ms = "./cmd_ms.inp"
     open(800,file=trim(ifile_cmd_ms),status="old")
     read(800,*) flag_ms_ff_LessPrint
     read(800,*) flag_fix_atoms_md
     read(800,*) flag_use_absorb_bound
     if(flag_use_absorb_bound) then
        backspace 800
        read(800,*) flag_use_absorb_bound,ABC_method !ABC_method="PML" or "Mur"
        if(ABC_method=="PML") then
           backspace 800
           read(800,*) flag_use_absorb_bound,ABC_method,R_pml,M_pml,nx_pml
           !usually, M_pml=2-4, R_pml=1d-4-1d-8, nx_pml=10-50
        endif
     endif
     read(800,*) flag_use_light_source
     if(flag_use_light_source) then
        backspace 800
        read(800,*) flag_use_light_source, ix_light_source
     endif

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

  call comm_bcast(flag_ms_ff_LessPrint, nproc_group_global)
  call comm_bcast(flag_fix_atoms_md,    nproc_group_global)
  call comm_bcast(flag_use_absorb_bound,nproc_group_global)
  call comm_bcast(flag_use_light_source,nproc_group_global)
  call comm_bcast(ix_light_source,      nproc_group_global)
  call comm_bcast(ABC_method,           nproc_group_global)
  call comm_bcast(R_pml,                nproc_group_global)
  call comm_bcast(M_pml,                nproc_group_global)
  call comm_bcast(nx_pml,               nproc_group_global)
  call comm_bcast(force_field_system,   nproc_group_global)
  call comm_bcast(nmol_s,               nproc_group_global)

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
  allocate( natom_mol(nmol)   )
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
    !read(800,*) dQave_thresh   !=1d-7 [au]  (from paper)
     read(800,*) conv_iter2     !=1d-10 [au] convergence of E(t) for pred-corrector
     read(800,*) a_ewald        != [1/A]
     read(800,*) r_cutoff       !=13.0 [A]   (from paper)
     read(800,*) G_cutoff       !=1.47 [1/A] (from paper)
     a_ewald  = a_ewald *au_length_aa !to [1/Bohr]
     r_cutoff = r_cutoff/au_length_aa !to [Bohr]
     G_cutoff = G_cutoff*au_length_aa !to [1/Bohr]
  endif

  call comm_bcast(iter_max,     nproc_group_global)
 !call comm_bcast(dQave_thresh, nproc_group_global)
  call comm_bcast(conv_iter2,   nproc_group_global)
  call comm_bcast(a_ewald,      nproc_group_global)
  call comm_bcast(r_cutoff,     nproc_group_global)
  call comm_bcast(G_cutoff,     nproc_group_global)

  !(read input no.5 & close the file)
  if (comm_is_root(nproc_id_global)) then
     read(800,*) tau_fric
!     read(800,*) flag_advanced  !=.false.
!     if(flag_advanced)then
!        read(800,*) keywd
!     endif
     close(800)
  endif

  call comm_bcast(tau_fric, nproc_group_global)

  !(temperature effect)
  tau_fric = tau_fric / au_time_fs  ![au] relax time constant
  if(tau_fric.le.1d-10) then
     fric = 0d0
     if(comm_is_root(nproc_id_global)) write(*,*) "#friction is not applied"
  else
     fric = 1d0/tau_fric
     if(comm_is_root(nproc_id_global)) write(*,*) "#friction is applied: fric=",real(fric)
  endif

  !(log)
  if (comm_is_root(nproc_id_global)) then
     if(flag_use_absorb_bound) then
        write(*,*) "flag_use_absorb_bound= ", flag_use_absorb_bound
        write(*,*) "ABC_method= ", trim(ABC_method)
        if(ABC_method=="PML") then
           write(*,*) "R, M, nx=", real(R_pml), M_pml, nx_pml
        endif
     endif
     if(flag_fix_atoms_md) write(*,*) "Atom is fixed"
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

  !allocation for multi-scale
  allocate( Qai_ms(4,nmol,nmacro_s:nmacro_e) )
  allocate( Qai_old_ms(4,nmol,nmacro_s:nmacro_e) )
  allocate( Qai_bk_ms(4,nmol,nmacro_s:nmacro_e) )
  allocate( Rxyz_wX_ms(3,4,nmol,nmacro_s:nmacro_e) )
  allocate( Rxyz_wX_ms_mid(3,4,nmol,nmacro_s:nmacro_e) )
  allocate( Rxyz_wX_ms_bk( 3,4,nmol,nmacro_s:nmacro_e) )
  allocate( Uene_ms(nmacro_s:nmacro_e) )
  allocate( Tene_ms(nmacro_s:nmacro_e) )
  allocate( Eall_ms(nmacro_s:nmacro_e) )
  allocate( Uene_intra_ms(nmacro_s:nmacro_e) )
  allocate( force_intra_m(3,NI,nmacro_s:nmacro_e) )
  if(flag_method_newton_Et) then
     allocate( vpr( 3,NI_wX,nmacro_s:nmacro_e), vmr(3,NI_wX,nmacro_s:nmacro_e) )
     allocate( rvec(3,NI_wX,nmacro_s:nmacro_e) )
     allocate( Qvec_bk(NI_wX,nmacro_s:nmacro_e) )
     allocate( vec2_tmp_m(NI_wX,nmacro_s:nmacro_e) )
     allocate( vec3_tmp_m(NI_wX,nmacro_s:nmacro_e) )
  endif


  !read coordinate of equiblium position (necessary?? maybe not)
  allocate( Rion_eq0(3,NI) )
  if (comm_is_root(nproc_id_global)) then     
     select case(iflag_atom_coor)
     case(ntype_atom_coor_cartesian)
        open(801,file='.atomic_coor.tmp',status='old')
        if(unit_length=='AA')then
           uconv = au_length_aa
        else  !au
           uconv = 1d0
        endif
        do ia = 1,NI
           read(801,*) char_atom, (tmpr(j),j=1,3),ik
           Rion_eq0(:,ia) = tmpr(:)/uconv
        enddo
        
     case(ntype_atom_coor_reduced)
        open(801,file='.atomic_red_coor.tmp',status='old')
        do ia = 1,NI
           read(801,*) char_atom, (tmpr(j),j=1,3),ik
           Rion_eq0(1,ia) = tmpr(1) * aLx
           Rion_eq0(2,ia) = tmpr(2) * aLy
           Rion_eq0(3,ia) = tmpr(3) * aLz
        enddo
        
     end select
         
     close(801)
  endif

  call comm_bcast(Rion_eq0,nproc_group_global)

  return
  
end subroutine init_cmd_maxwell_ms

subroutine cmd_maxwell_ms
  use Global_Variables
  use timer
  use opt_variables
  use performance_analyzer
  use salmon_parallel
  use salmon_communication
  use salmon_file
  use misc_routines
  use inputoutput, only: t_unit_time, t_unit_current, t_unit_elec
  implicit none
  logical :: flag_iter2_loop
  integer :: iter, ix_m, iy_m, iz_m, imacro,igroup,i, index
  logical :: flg_out_ms_step, flg_out_ms_next_step
  real(8) :: Temperature_ion, gkT, Qnh
  real(8),allocatable :: Enh_m(:), Enh_gkTlns_m(:)
  real(8),allocatable :: Mt_ms(:,:)
  character(100) :: comment_line

  iter2=1

#ifdef ARTED_LBLK
  call opt_vars_init_t4ppt()
#endif

  !(unnecessary)
  deallocate( dRion, dRion_m )
  deallocate( javt )

  !
  allocate( Mt_ms(3,nmacro_s:nmacro_e) )
  allocate( Enh_m(nmacro_s:nmacro_e) )
  allocate( Enh_gkTlns_m(nmacro_s:nmacro_e) )

  !if material region does not exist, nmacro-->0, nmacro_e-->0
  if(nx_m .le. 0)then
     nmacro   = 0
     nmacro_e = 0
     if(comm_is_root(nproc_id_global)) then
        write(*,*)"# no material region: set nmacro-->0, nmacro_e-->0"
     endif
  endif

  do imacro = nmacro_s, nmacro_e
     call take_back_mol_into_ucell_ms(imacro)
  enddo


  if(comm_is_root(nproc_id_global)) then
     write(*,*) 'This is the end of preparation for Real time calculation'
     call timer_show_current_hour('elapse time=',LOG_ALL)
     write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation============================

  if(flag_use_light_source) then
     if(comm_is_root(nproc_id_global)) write(*,*) "use light source"
  endif

  if(flag_use_absorb_bound .and. ABC_method=="PML") &
       call init_absorb_bound_1d_PML
    
  deallocate(zu_t)

  !(for restarting with read_rt_wfn_k_ms=y option)
  if(read_rt_wfn_k_ms=='y') then
     call add_init_Ac_ms
     call assign_mp_variables_omp_EB()
     do imacro = 1, nmacro
        ix_m = macropoint(1,imacro)
        iy_m = macropoint(2,imacro)
        iz_m = macropoint(3,imacro)
        Jm_m(1:3,imacro)=matmul(trans_mat(1:3,1:3),Jm_ms(1:3,ix_m,iy_m,iz_m))
        Jm_new_m(1:3,imacro)=matmul(trans_mat(1:3,1:3),Jm_new_ms(1:3,ix_m,iy_m,iz_m))
     end do
  endif

!reentrance

  !position_option='rewind'
  entrance_iter=-1
  call reset_rt_timer

  ! create directory for exporting
  call create_dir_ms

  do imacro = nmacro_s, nmacro_e
     call cal_force_energy_CRK_MS(imacro)
     if(read_rt_wfn_k_ms=='n') call init_some_var_ms(imacro)
  enddo

  ! Export to file_trj (initial step)
  if (out_rvf_rt=='y')then
      write(comment_line,110) -1, 0.0d0

      do imacro = nmacro_s, nmacro_e
         if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
            & write(comment_line,112) trim(comment_line), xi_nh_m(imacro)

         if(flag_ms_ff_LessPrint)then
            if (imacro==1)then
            call write_xyz_ms(comment_line,"new","rvf",imacro)
            call write_xyz_ms(comment_line,"add","rvf",imacro)
            endif
         else
            call write_xyz_ms(comment_line,"new","rvf",imacro)
            call write_xyz_ms(comment_line,"add","rvf",imacro)
         endif
      enddo
      call comm_sync_all
  endif

  !(get initial ion current for initial step)
  if(read_rt_wfn_k_ms=='n') then
     do imacro = nmacro_s, nmacro_e
        Mt_ms(:,imacro) =0d0  !for initial step
        call cal_dipole_moment_ms(Mt_ms(:,imacro),dt,imacro)
        jav(:) = 0d0  !for initial step
        if (comm_is_root(nproc_id_tdks)) then
           jm_new_m_tmp(1:3,imacro) = jav(1:3)
        end if
     enddo
     call comm_summation(jm_new_m_tmp,Jm_m,3*nmacro,nproc_group_global)
     do imacro = 1, nmacro
        ix_m = macropoint(1,imacro)
        iy_m = macropoint(2,imacro)
        iz_m = macropoint(3,imacro)
        Jm_ms(1:3,ix_m,iy_m,iz_m) = matmul(trans_inv(1:3,1:3),Jm_m(1:3,imacro))
     enddo
  endif

  Enh_gkTlns_m(:) = 0d0
  Enh_m(:)        = 0d0
  if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
     gkT = 3d0*NI * kB/hartree2J*temperature0_ion
     Qnh = gkT * thermostat_tau**2d0
  endif


  RTiteratopm : do iter=entrance_iter+1, Nt

    !! NOTE: flg_out_ms_step (the macroscopic field will exported in this step)
    flg_out_ms_step = .false.
    if (mod(iter, out_ms_step)==0) flg_out_ms_step=.true.

    iter_save=iter

    !! Update of the Microscopic System
    flg_out_ms_next_step=.false.
    if( mod(iter+1, out_ms_step)==0 ) flg_out_ms_next_step=.true.

    iter2=1
    flag_conv_iter2(:)=.false.

10  continue

    Macro_loop : do imacro = nmacro_s, nmacro_e
      
      if(iter2==1) then
         !! update coordinate and velocity in MD option (part-1)
         Rion_m_mid(:,:,imacro) = Rion_m(:,:,imacro)
         if(.not.flag_fix_atoms_md) &
         call dt_evolve_MD_1_MS(iter,imacro)
         velocity_m_mid(:,:,imacro) = velocity_m(:,:,imacro)
         Rion_m_mid(:,:,imacro)= 0.5d0*(Rion_m_mid(:,:,imacro) + Rion_m(:,:,imacro))
      else
         velocity_m(:,:,imacro)=velocity_m_mid(:,:,imacro)
      endif !iter2

      ! force
      if(iter2==1) then
         Rxyz_wX_ms_bk(:,:,:,imacro) = Rxyz_wX_ms(:,:,:,imacro) 
         Qai_bk_ms(:,:,imacro) = Qai_ms(:,:,imacro)
      endif
      call cal_force_energy_CRK_MS(imacro)
      Rxyz_wX_ms_mid(:,:,:,imacro) = 0.5d0*(Rxyz_wX_ms_bk(:,:,:,imacro)+Rxyz_wX_ms(:,:,:,imacro))

      !! update coordinate and velocity in MD option (part-2)
      if(.not.flag_fix_atoms_md) &
      call dt_evolve_MD_2_MS(iter,imacro)

      ! current
      call cal_dipole_moment_current_ms_EB(Mt_ms(:,imacro),jav,dt,imacro)

      if(flag_method_newton_Et) call save_matrix_for_newton_method(imacro)

      if(Sym/=1) jav(1:2)=0d0
      !! Special rule for debug mode (macroscopic system does not have a matter)
      !! NOTE: This condition will be removed and replaced by FDTD mode
      !!       after merging the common Maxwell calculation routine.
      if (debug_switch_no_radiation) jav(:)=0d0

      jav(1)=0d0  !force to zero for propagation direction
      if (comm_is_root(nproc_id_tdks)) then
         jm_new_m_tmp(1:3,imacro) = jav(1:3)
      end if

      ! Calculate + store electron energy
      if (flg_out_ms_next_step) then !! mod(iter+1,out_ms_step) == 0
        if(comm_is_root(nproc_id_tdks)) then
          energy_elec_Matter_new_m_tmp(imacro) = Eall_ms(imacro) - Eall0_m(imacro)
        end if
      end if

    end do Macro_loop !end of Macro_loop iteraction========================

    !===========================================================================

    call comm_summation(jm_new_m_tmp, Jm_new_m, 3 * nmacro, nproc_group_global)

    if(flg_out_ms_next_step) then
       call comm_summation(energy_elec_Matter_new_m_tmp, energy_elec_Matter_new_m, nmacro, nproc_group_global)
    end if
    
!$omp parallel do default(shared) private(imacro,ix_m,iy_m,iz_m)
    do imacro = 1, nmacro
      ix_m = macropoint(1,imacro)
      iy_m = macropoint(2,imacro)
      iz_m = macropoint(3,imacro)
      Jm_new_ms(1:3,ix_m,iy_m,iz_m) = matmul(trans_inv(1:3,1:3), Jm_new_m(1:3,imacro))
    end do
!$omp end parallel do


    call proceed_ms_variables_omp_EB(iter2)
    call dt_evolve_EB_1d_cmd() 
    call assign_mp_variables_omp_EB()

    iter2=iter2+1
    if(iter2==100) then
       write(*,*) "iteration of predictor-corrector did not converge",iter2
       stop
    endif

    flag_iter2_loop = .false.
    do imacro=1,nmacro
       dEt_pc = sqrt(sum((E_m_bk(:,imacro) - E_m(:,imacro))**2))
       if( dEt_pc .lt. conv_iter2) then
          flag_conv_iter2(imacro)=.true.  !not used now
       else
          flag_conv_iter2(imacro)=.false.
          flag_iter2_loop=.true.
       endif
    enddo

    call comm_sync_all
    if(flag_iter2_loop) goto 10

    !! Compute EM field energy and other important quantities
    call calc_energy_joule_EB()
    if(flag_use_light_source) &
    call calc_energy_joule_subtract_light_source_EB() !subtract light source J
    if (flg_out_ms_step) then !! mod(iter, out_ms_step) == 0
       call calc_energy_elemag_EB()
       call calc_total_energy_EB
    end if

    !! Store data_vac_ac variables
    call store_data_local_ac_jm()
    call store_data_vac_ac()
    !! Store data_out variables
    if (flg_out_ms_step) then !! mod(iter, out_ms_step) == 0
       index = iter / out_ms_step
       if(.not. flag_ms_ff_LessPrint) then !AY
          call store_data_out_omp(index)
          call write_data_out(index) !now moved to here
       endif
    end if

    !! Print log (energy)
    if(flg_out_ms_step .and. comm_is_root(nproc_id_global)) then
       call trace_ms_calculation()
    end if


    !===========================================================================
    ! Export to file_trj
    if (out_rvf_rt=='y' .and. mod(iter,out_rvf_rt_step)==0)then
       write(comment_line,110) iter, iter*dt
110    format("#rt   step=",i8,"   time",e16.6)
112    format(a,"  xi_nh=",e18.10)
       if(flag_ms_ff_LessPrint) then  !AY
          do imacro = nmacro_s, nmacro_e
             if(imacro==1) then
                call write_xyz_ms(comment_line,"add","rvf",imacro)
                if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
                &  write(comment_line,112) trim(comment_line),xi_nh_m(imacro)
             endif
          enddo
          call comm_sync_all
       else   !AY
          do igroup=1,ndivide_macro
             do i=1,nmacro_write_group
                imacro = (igroup-1)*nmacro_write_group + i
                if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e) then
                   if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
                   &  write(comment_line,112) trim(comment_line),xi_nh_m(imacro)
                   call write_xyz_ms(comment_line,"add","rvf",imacro)
                endif
             enddo
             call comm_sync_all
          enddo
       endif
    endif

    !===========================================================================
    
  enddo RTiteratopm !end of RT iteraction========================


  if(comm_is_root(nproc_id_global)) then
    write(*,'(1X,A)') '-----------------------------------------------'
  end if
  call write_performance(trim(directory)//'ms_performance')

  if(comm_is_root(nproc_id_global)) write(*,*) 'This is the start of write section'

  !!!xxxx move into the do loop of time later hoge 
  !! Write data out by using MPI
  !if(.not. flag_ms_ff_LessPrint) then !AY
  !   do index = 0, Ndata_out-1
  !      call write_data_out(index)
  !   end do
  !endif  
  call write_data_local_ac_jm()
  call write_data_vac_ac()

  ! Export last atomic coordinate and velocity & Close file_trj
  if(flag_ms_ff_LessPrint) then !AY
     do imacro = nmacro_s, nmacro_e
        if(imacro==1) then
           call write_xyz_ms(comment_line,"end","rvf",imacro)
           call write_ini_coor_vel_ms_for_restart(imacro)
        endif
     enddo
     call comm_sync_all
  else  !AY
     do igroup=1,ndivide_macro
        do i=1,nmacro_write_group
           imacro = (igroup-1)*nmacro_write_group + i
           if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e) then
              call write_xyz_ms(comment_line,"end","rvf",imacro)
              call write_ini_coor_vel_ms_for_restart(imacro)
           endif
        enddo
        call comm_sync_all
     enddo
  endif


  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of write section'
  end if

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of RT calculation'
    write(*,*) '-----------------------------------------------------------'
  end if

  call comm_sync_all

  if(comm_is_root(nproc_id_global)) write(*,*) 'This is the end of all calculation'
1 if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'

!========================================================
contains

  subroutine reset_rt_timer
    implicit none
    integer :: i
    do i = LOG_DT_EVOLVE,LOG_ALLREDUCE
      call timer_reset(i)
    end do
  end subroutine
  
  subroutine proceed_ms_variables_omp_EB(iflag)
    implicit none
    integer :: iflag, iimacro, iix_m, iiy_m, iiz_m, ipt

    iy_m = ny_origin_m
    iz_m = nz_origin_m

    if(iflag==1) then

       do ipt = npt_intpl_E-1,2,-1
          E_hist_m(:,:,ipt) = E_hist_m(:,:,ipt-1)
       enddo
       E_hist_m(:,:,1) = E_m(:,:)
       if(flag_use_absorb_bound)then
       if(ABC_method=="Mur") then
          E_last(:,nx1_m:nx1_m+1,iy_m,iz_m) = E_ms(:,nx1_m:nx1_m+1,iy_m,iz_m)
          E_last(:,nx2_m-1:nx2_m,iy_m,iz_m) = E_ms(:,nx2_m-1:nx2_m,iy_m,iz_m)
          B_last(:,nx1_m:nx1_m+1,iy_m,iz_m) = B_ms(:,nx1_m:nx1_m+1,iy_m,iz_m)
          B_last(:,nx2_m-1:nx2_m,iy_m,iz_m) = B_ms(:,nx2_m-1:nx2_m,iy_m,iz_m)
       else if(ABC_method=="PML") then
          E_last(:,:,iy_m,iz_m) = E_ms(:,:,iy_m,iz_m)
          B_last(:,:,iy_m,iz_m) = B_ms(:,:,iy_m,iz_m)
       endif
       endif

    else
       E_ms(:,:,:,:) = E_last(:,:,:,:) 
       B_ms(:,:,:,:) = B_last(:,:,:,:)
    endif

!$omp parallel do collapse(3) default(shared) private(iix_m, iiy_m, iiz_m)
    do iiz_m = nz1_m, nz2_m
    do iiy_m = ny1_m, ny2_m
    do iix_m = nx1_m, nx2_m
       Jm_old_ms(1:3,iix_m,iiy_m,iiz_m) = Jm_ms    (1:3,iix_m,iiy_m,iiz_m)
       Jm_ms    (1:3,iix_m,iiy_m,iiz_m) = Jm_new_ms(1:3,iix_m,iiy_m,iiz_m)
       Jm_new_ms(1:3,iix_m,iiy_m,iiz_m) = 0d0
    end do
    end do
    end do
!$omp end parallel do

!$omp parallel do default(shared) private(iimacro, iix_m, iiy_m, iiz_m)
    do iimacro = 1, nmacro
      Jm_m(1:3, iimacro) = Jm_new_m(1:3, iimacro)
      if (flg_out_ms_step) then
        iix_m = macropoint(1, iimacro)
        iiy_m = macropoint(2, iimacro)
        iiz_m = macropoint(3, iimacro)
        energy_elec_ms(iix_m, iiy_m, iiz_m) = energy_elec_Matter_new_m(iimacro)
      end if
    end do
!$omp end parallel do

  end subroutine proceed_ms_variables_omp_EB

  subroutine assign_mp_variables_omp_EB()
    implicit none
    integer :: iix_m, iiy_m, iiz_m, iimacro

!$omp parallel do default(shared) private(iimacro,iix_m,iiy_m,iiz_m)
    do iimacro = 1, nmacro
       iix_m = macropoint(1,iimacro)
       iiy_m = macropoint(2,iimacro)
       iiz_m = macropoint(3,iimacro)
       !! Assign the E and H field into the local macropoint variables
       E_m_bk(1:3,iimacro) = E_m(1:3,iimacro)
       E_m(1:3,iimacro) = matmul(trans_mat(1:3,1:3),E_ms(1:3,iix_m,iiy_m,iiz_m))
       B_m(1:3,iimacro) = matmul(trans_mat(1:3,1:3),B_ms(1:3,iix_m,iiy_m,iiz_m))
    end do
!$omp end parallel do

  end subroutine assign_mp_variables_omp_EB
  
  
  subroutine store_data_local_ac_jm()
    implicit none
    integer :: iimacro
!$omp parallel do default(shared) private(iimacro)
    do iimacro = nmacro_s, nmacro_e
      !! Store data_local_Ac, data_local_Jm
      data_local_Ac(1:3,iimacro,iter) = E_m(1:3,iimacro)
      data_local_jm(1:3, iimacro, iter) = Jm_m(1:3, iimacro)
    end do
!$omp end parallel do
  end subroutine store_data_local_ac_jm
  

  subroutine store_data_out_omp(index)
    use Global_Variables
    implicit none
    integer, intent(in) :: index
    integer :: iproc, ipos
    integer :: iix_m, iiy_m, iiz_m
    
    iproc = mod(index, nproc_size_global)
    ipos = (index - iproc) / nproc_size_global
    if (iproc == nproc_id_global) then
!$omp parallel do collapse(3) default(shared) private(iiy_m, iix_m, iiz_m)
        do iiz_m = nz1_m, nz2_m
        do iiy_m = ny1_m, ny2_m
        do iix_m = nx1_m, nx2_m
            data_out(1:3,  iix_m,iiy_m,iiz_m,ipos) = Ac_ms(1:3,iix_m,iiy_m,iiz_m)
            data_out(4:6,  iix_m,iiy_m,iiz_m,ipos) = E_ms(1:3, iix_m,iiy_m,iiz_m)
            data_out(7:9,  iix_m,iiy_m,iiz_m,ipos) = B_ms(1:3, iix_m,iiy_m,iiz_m)
            data_out(10:12,iix_m,iiy_m,iiz_m,ipos) =Jm_ms(1:3, iix_m,iiy_m,iiz_m)
            data_out(13,   iix_m,iiy_m,iiz_m,ipos) = Energy_elec_ms(iix_m,iiy_m,iiz_m)
            data_out(14,   iix_m,iiy_m,iiz_m,ipos) = Energy_joule_ms(iix_m,iiy_m,iiz_m)
            data_out(15,   iix_m,iiy_m,iiz_m,ipos) = Energy_elemag_ms(iix_m,iiy_m,iiz_m)
        end do
        end do
        end do
!$omp end parallel do
    end if
  end subroutine
  
  
  subroutine create_dir_ms
    implicit none
    integer :: imacro

    call create_directory(dir_ms)
    call create_directory(dir_ms_RT_Ac)

    do imacro = nmacro_s, nmacro_e
       call create_directory(dir_ms_M(imacro))
    enddo

  end subroutine

  subroutine take_back_mol_into_ucell_ms(imacro)
    implicit none
    integer j,imol,iatom_c,iatom,ia, imacro

    do imol=1,nmol
       iatom_c = mol2atom_cnt(imol)
       do j=1,3
          if( Rion_m(j,iatom_c,imacro).lt.0d0 ) then
             do ia= 1,natom_mol(imol)
                iatom = mol2atom_top(imol) + ia-1
                Rion_m(j,iatom,imacro) = Rion_m(j,iatom,imacro) + aL(j)
             enddo
          else if( Rion_m(j,iatom_c,imacro).gt.aL(j) ) then
             do ia= 1,natom_mol(imol)
                iatom = mol2atom_top(imol) + ia-1
                Rion_m(j,iatom,imacro) = Rion_m(j,iatom,imacro) - aL(j)
             enddo
          endif
       enddo
    enddo

  end subroutine

  subroutine write_data_out(index)
    implicit none
    integer, intent(in) :: index
    integer :: iix_m, iiy_m, iiz_m, fh_ac
    integer :: iproc, ipos
    
    iproc = mod(index, nproc_size_global)
    ipos  = (index - iproc) / nproc_size_global

    if (iproc == nproc_id_global) then
      write(file_ac, "(A,A,'_Ac_',I6.6,'.data')") & 
        & trim(dir_ms_RT_Ac), trim(SYSname),  out_ms_step * index
      fh_ac = open_filehandle(file_ac)

      write(fh_ac, '("#",1X,A)') "Macroscopic field distribution"
      write(fh_ac, '("#",1X,A,":",1X,A)') "IX,IY,IZ", "Coordinate"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Ac", "Vector potential field"
      write(fh_ac, '("#",1X,A,":",1X,A)') "E", "Electric field"
      write(fh_ac, '("#",1X,A,":",1X,A)') "B", "Magnetic field"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Eex", "Electron excitation energy"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Eabs", "Absorbed energy"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Eemf", "Total EM field energy"

      write(fh_ac, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 1,  "Ix", "none", &
        & 2,  "Iy", "none", &
        & 3,  "Iz", "none", &
        & 4,  "Ac_x", "a.u.", & !!, trim(t_unit_ac%name), &
        & 5,  "Ac_y", "a.u.", & !!, trim(t_unit_ac%name), &
        & 6,  "Ac_z", "a.u.", & !!, trim(t_unit_ac%name), &
        & 7,  "E_x",  "a.u.", & !!, trim(t_unit_current%name), &
        & 8,  "E_y",  "a.u.", & !!, trim(t_unit_current%name), &
        & 9,  "E_z",  "a.u.", & !!, trim(t_unit_current%name), &
        & 10, "B_x",  "a.u.", & !!, trim(t_unit_current%name), &
        & 11, "B_y",  "a.u.", & !!, trim(t_unit_current%name), &
        & 12, "B_z",  "a.u.", & !!, trim(t_unit_current%name), &
        & 13, "Jm_x", "a.u.", & !!, trim(t_unit_current%name), &
        & 14, "Jm_y", "a.u.", & !!, trim(t_unit_current%name), &
        & 15, "Jm_z", "a.u.", & !!, trim(t_unit_current%name), &
        & 16, "Eex",  "a.u./unitcell", & !!, trim(t_unit_current%name), &
        & 17, "Eabs", "a.u./unitcell", & !!, trim(t_unit_current%name), &
        & 18, "Eemf", "a.u./unitcell" !!, trim(t_unit_current%name)
      write(fh_ac,*)

      !! TODO: Support the automatic unit-system conversion of _ac.data files
      do iiz_m = nz1_m, nz2_m
      do iiy_m = ny1_m, ny2_m
      do iix_m = nx1_m, nx2_m
         write(fh_ac,'(I6,1X,I6,1X,I6,99(1X,E23.15E3))',advance='no')  &
              & iix_m, iiy_m, iiz_m, &
              & data_out(1:ndata_out_column, iix_m, iiy_m, iiz_m, ipos)
         write(fh_ac,*)
      end do
      end do
      end do
      close(fh_ac)
    end if
   
    return
  end subroutine write_data_out


  subroutine store_data_vac_ac()
    implicit none
    ! Export the Ac field of detecting point
    if(comm_is_root(nproc_id_global)) then
       data_vac_Ac(1:3, 1, iter) = E_ms(1:3,ix_detect_l,iy_detect,iz_detect)
       data_vac_Ac(1:3, 2, iter) = E_ms(1:3,ix_detect_r,iy_detect,iz_detect)
    end if
  end subroutine store_data_vac_ac
  
  
  subroutine trace_ms_calculation()
    implicit none
    integer :: iix_m, iiy_m, iiz_m
    real(8) :: rrx, rry, rrz, eem, sem
    real(8) :: sx1, sy1, sz1
    real(8) :: sx2, sy2, sz2
    real(8) :: cx1, cy1, cz1
    real(8) :: cx2, cy2, cz2
    real(8) :: wx, wy, wz
    
    write(*,'(1X,A)') '-----------------------------------------------'
    write(*,'(1X,A,I6,A,I6)') 'Multiscale iter =', iter, '/', Nt
    
    write(*,'(1X,A)') 'Microscopic system:'
    
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Excitation energy:', &
    & "Eex =", total_energy_elec, &
    & "Diff =", total_energy_elec - total_energy_elec_old, &
    & "[au]"

    write(*,'(1X,A)') 'Macroscopic system:'
    
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Absorbed energy:', &
      & "Eabs =", total_energy_absorb, &
      & "Diff =", total_energy_absorb - total_energy_absorb_old, &
      & "[au]"
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Field energy:', &
      & "Eemf =", total_energy_elemag, &
      & "Diff =", total_energy_elemag - total_energy_elemag_old, &
      & "[au]"
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Total EM energy:', &
      & "Etot =", total_energy_em, &
      & "Diff =", total_energy_em - total_energy_em_old, &
      & "[au]"
      
    sem = 0d0;
    sx1 = 0d0; sy1 = 0d0; sz1 = 0d0
    sx2 = 0d0; sy2 = 0d0; sz2 = 0d0
!$omp parallel do default(shared) collapse(3) &
!$omp private(iix_m, iiy_m, iiz_m, rrx, rry, rrz, eem) &
!$omp reduction(+: sx1, sy1, sz1, sx2, sy2, sz2, sem)
    do iiz_m = nz1_m, nz2_m
      do iiy_m = ny1_m, ny2_m
        do iix_m = nx1_m, nx2_m
          rrx = iix_m * HX_m
          rry = iiy_m * HY_m
          rrz = iiz_m * HZ_m
          
          eem = Energy_elemag_ms(iix_m, iiy_m, iiz_m)
          sem = sem + eem
          
          sx1 = sx1 + eem * rrx
          sy1 = sy1 + eem * rry
          sz1 = sz1 + eem * rrz
          
          sx2 = sx2 + eem * rrx ** 2
          sy2 = sy2 + eem * rry ** 2
          sz2 = sz2 + eem * rrz ** 2
        end do
      end do
    end do
!$omp end parallel do

    if (0d0 < sem) then
      cx1 = sx1 / sem; cy1 = sy1 / sem; cz1 = sz1 / sem 
      cx2 = sx2 / sem; cy2 = sy2 / sem; cz2 = sz2 / sem 
    else
      cx1 = 0d0; cy1 = 0d0; cz1 = 0d0
      cx2 = 0d0; cy2 = 0d0; cz2 = 0d0
    endif
    wx = sqrt(abs(cx2-cx1**2))
    wy = sqrt(abs(cy2-cy1**2))
    wz = sqrt(abs(cz2-cz1**2))
    write(*, '(1X,A)') "Position of wavepacket:"
    write(*,'(2X,A20,A10,3(ES12.3E3,","),A5)') &
      & 'Central position:', "<r> =", cx1, cy2, cz2, "[au]"
    write(*,'(2X,A20,A10,3(ES12.3E3,","),A5)') &
      & 'Spatial spreading:', "rms =", wx, wy, wz, "[au]"
    return 
  end subroutine trace_ms_calculation


  !===============================================================

  subroutine write_data_local_ac_jm()
    use salmon_file
    implicit none
    integer :: fh_ac_m
    integer :: iimacro, iiter
    character(100) :: file_ac_m

    !do iimacro = nmacro_s, nmacro_e
    if(flag_ms_ff_LessPrint) then !AY

       do iimacro = nmacro_s, nmacro_e

          if(iimacro.ne.1) cycle

          if(comm_is_root(nproc_id_tdks)) then

          if(iimacro.ge.nmacro_s .and. iimacro.le.nmacro_e) then

          write(file_ac_m,"(A,A,'_Ac_M.data')")trim(dir_ms_M(iimacro)),trim(SYSname)
          fh_ac_m = open_filehandle(file_ac_m)

          write(fh_ac_m, '("#",1X,A)') "Local variable at macro point"
        
          write(fh_ac_m, "('#',1X,A,':',3(1X,I6))") "Macropoint", macropoint(1:3, iimacro)
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "E", "External electric field"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Tmp", "Temperature"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Uene", "Potential energy(including field interaction)/(unit cell)"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Tene", "Kinetic energy/(unit cell)"

          write(fh_ac_m, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
                  & 1, "Time", trim(t_unit_time%name), &
                  & 2, "E_x", trim(t_unit_elec%name), &
                  & 3, "E_y", trim(t_unit_elec%name), &
                  & 4, "E_z", trim(t_unit_elec%name), &
                  & 5, "Jm_x", trim(t_unit_current%name), &
                  & 6, "Jm_y", trim(t_unit_current%name), &
                  & 7, "Jm_z", trim(t_unit_current%name), &
                  & 8, "Tmp",       "K",  &
                  & 9, "Uene",      "au", &
                  & 10,"Uene+Tene", "au"
          write(fh_ac_m,*)
          do iiter = 0, Nt
             write(fh_ac_m, "(F16.8,9(1X,ES22.14E3,1X))",advance='no') &
                     & iiter * dt * t_unit_time%conv, &
                     & data_local_Ac(1:3, iimacro, iiter) * t_unit_elec%conv, &
                     & data_local_jm(1:3, iimacro, iiter) * t_unit_current%conv, &
                     & data_local_MD(1:3, iimacro, iiter)
             write(fh_ac_m,*)
          end do !iiter

          close(fh_ac_m)

          endif !<--if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e)
          endif ! nproc_id_tdks

       enddo !iimacro

       call comm_sync_all

    else

    do igroup=1,ndivide_macro
       do i=1,nmacro_write_group

          if(comm_is_root(nproc_id_tdks)) then

          iimacro = (igroup-1)*nmacro_write_group + i
          if(iimacro.ge.nmacro_s .and. iimacro.le.nmacro_e) then

          write(file_ac_m,"(A,A,'_Ac_M.data')")trim(dir_ms_M(iimacro)),trim(SYSname)
          fh_ac_m = open_filehandle(file_ac_m)

          write(fh_ac_m, '("#",1X,A)') "Local variable at macro point"
          write(fh_ac_m, "('#',1X,A,':',3(1X,I6))") "Macropoint", macropoint(1:3, iimacro)

          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "E", "External electric field"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Tmp", "Temperature"

          write(fh_ac_m, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
                  & 1, "Time", trim(t_unit_time%name), &
                  & 2, "E_x", trim(t_unit_elec%name), &
                  & 3, "E_y", trim(t_unit_elec%name), &
                  & 4, "E_z", trim(t_unit_elec%name), &
                  & 5, "Jm_x", trim(t_unit_current%name), &
                  & 6, "Jm_y", trim(t_unit_current%name), &
                  & 7, "Jm_z", trim(t_unit_current%name), &
                  & 8, "Tmp",       "K",  &
                  & 9, "Uene",      "au", &
                  & 10,"Uene+Tene", "au"
          write(fh_ac_m,*)
          do iiter = 0, Nt
             write(fh_ac_m, "(F16.8,9(1X,E23.15E3,1X))",advance='no') &
                     & iiter * dt * t_unit_time%conv, &
                     & data_local_Ac(1:3, iimacro, iiter) * t_unit_elec%conv, &
                     & data_local_jm(1:3, iimacro, iiter) * t_unit_current%conv, &
                     & data_local_MD(1:3, iimacro, iiter)
             write(fh_ac_m,*)
          end do !iiter

          close(fh_ac_m)

          endif !<--if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e)
          endif ! nproc_id_tdks

       enddo !i
       call comm_sync_all
    enddo !igroup
    endif

  end subroutine
  
  
  subroutine write_data_vac_ac()
    use salmon_file
    implicit none
    integer :: fh_ac_vac
    integer :: iiter
    character(100) :: file_ac_vac
    
    if (comm_is_root(nproc_id_global)) then
      write(file_ac_vac,"(A, A, '_Ac_vac.data')") trim(dir_ms_RT_Ac),trim(SYSname)
      fh_ac_vac = open_filehandle(file_ac_vac)

      write(fh_ac_vac, '("#",1X,A)') "E(t) vacuum region"
      write(fh_ac_vac, '("#",1X,A)') "Data of E(t) field at the end of media"
      
      write(fh_ac_vac, "('#',1X,A,':',3(1X,I6))") "L", ix_detect_l, iy_detect, iz_detect
      write(fh_ac_vac, "('#',1X,A,':',3(1X,I6))") "R", ix_detect_r, iy_detect, iz_detect
      
      write(fh_ac_vac, '("#",99(1X,I0,":",A,"[",A,"]"))') &
              & 1, "Time", trim(t_unit_time%name), &
              & 2, "E_x(L)", trim(t_unit_elec%name), &
              & 3, "E_y(L)", trim(t_unit_elec%name), &
              & 4, "E_z(L)", trim(t_unit_elec%name), &
              & 5, "E_x(R)", trim(t_unit_elec%name), &
              & 6, "E_y(R)", trim(t_unit_elec%name), &
              & 7, "E_z(R)", trim(t_unit_elec%name)
      do iiter = 0, Nt
         write(fh_ac_vac, "(F16.8,6(1X,E23.15E3,1X))") &
                 & iiter * dt * t_unit_time%conv, &
                 & data_vac_Ac(1:3, 1, iiter) * t_unit_elec%conv, &
                 & data_vac_Ac(1:3, 2, iiter) * t_unit_elec%conv
      end do

      close(fh_ac_vac)
    end if
    call comm_sync_all  
  end subroutine write_data_vac_ac

  subroutine write_xyz_ms(comment,action,rvf,imacro)

  ! Write xyz in xyz format but also velocity and force are printed if necessary
  ! (these can be used for restart of opt and md)
    use Global_Variables
    use inputoutput, only: au_length_aa
    use salmon_global, only: SYSname,iflag_atom_coor,ntype_atom_coor_cartesian,ntype_atom_coor_reduced
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_file
    implicit none
    integer :: ia,unit_xyz,imacro 
    character(3) :: action,rvf
    character(1024) :: file_trj
    character(*) :: comment

    if(comm_is_root(nproc_id_tdks)) then

      unit_xyz = 1000 + imacro

      if(action=='new') then

        write(file_trj, "(A,A,'_trj.xyz')") trim(dir_ms_M(imacro)),trim(SYSname)
        open(unit_xyz,file=trim(file_trj),status="unknown")

      else if(action=='add') then

         write(unit_xyz,*) NI
         write(unit_xyz,*) trim(comment)
         do ia=1,NI
            if(      rvf=="r  " ) then
               write(unit_xyz,100) trim(atom_name(ia)), &
                           & Rion_m(1:3,ia,imacro)*au_length_aa
            else if( rvf=="rv " ) then
               write(unit_xyz,110) trim(atom_name(ia)), &
                           & Rion_m(1:3,ia,imacro)*au_length_aa,&
                           & velocity_m(1:3,ia,imacro)
            else if( rvf=="rvf" ) then
               write(unit_xyz,120) trim(atom_name(ia)), &
                           & Rion_m(1:3,ia,imacro)*au_length_aa,&
                           & velocity_m(1:3,ia,imacro), force_m(1:3,ia,imacro)
            endif
         enddo

      else if(action=='end') then

         close(unit_xyz)

      endif

100 format(a2,3f18.10)
110 format(a2,3f18.10, "  #v=",3f18.10)
120 format(a2,3f18.10, "  #v=",3f18.10, "  #f=",3f18.10)

    endif

!    call comm_sync_all
  end subroutine

  subroutine write_ini_coor_vel_ms_for_restart(imacro)
  !print atomic coordinate and velocity at the last step just to set_ini_coor_vel option
    use Global_Variables
    use inputoutput, only: au_length_aa
!    use salmon_global, only: SYSname,iflag_atom_coor,ntype_atom_coor_cartesian,ntype_atom_coor_reduced
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_file
    implicit none
    integer :: ia,unit_ini_rv,imacro 
    character(1024) :: file_ini_rv
    real(8) :: ulength_to_au

    if(comm_is_root(nproc_id_tdks)) then

      unit_ini_rv = 5000 + imacro

      write(file_ini_rv, "(A,'last_coor_vel.data')") trim(dir_ms_M(imacro))
      open(unit_ini_rv,file=trim(file_ini_rv),status="unknown")

      select case(unit_length)
      case('au','a.u.')
         ulength_to_au   = 1d0
      case('AA','angstrom','Angstrom')
         ulength_to_au   = 1d0/au_length_aa
      end select
      
      do ia=1,NI
         write(unit_ini_rv,110)  &
              & Rion_m(1:3,ia,imacro)/ulength_to_au,&
              & velocity_m(1:3,ia,imacro)
      enddo
110   format(6e24.12)

      close(unit_ini_rv)

    endif

!    call comm_sync_all
  end subroutine

  subroutine dt_evolve_MD_1_MS(iter,imacro)
  ! Velocity Verlet integrator for ion dynamics in multi-scale
    use md_ground_state
    implicit none
    integer :: iter,imacro,ia
    real(8) :: dt_h,mass_au
    character(100) :: comment_line

    dt_h = dt/2d0
    velocity(:,:) = velocity_m(:,:,imacro)
    force(:,:)    = force_m(:,:,imacro)
    Rion(:,:)     = Rion_m(:,:,imacro)
    Rion_eq(:,:)  = Rion_eq_m(:,:,imacro)

    !NHC act on velocity with dt/2
    if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
       call apply_nose_hoover_velocity_ms(imacro,dt_h)
    endif

    !update velocity with dt/2
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
      !velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h  &
       &                               - fric * velocity(:,ia) * dt_h
    enddo

    !velocity scaling
    if(step_velocity_scaling>=1 .and. mod(iter,step_velocity_scaling)==0)then
       call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity(:,:))
       call apply_velocity_scaling_ion(Temperature_ion,velocity(:,:))
    endif

    !update ion coordinate with dt
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
       Rion(:,ia) = Rion(:,ia) + velocity(:,ia)*dt
    enddo

    !NHC act on thermostat with dt
    if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
       call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
       call apply_nose_hoover_thermostat_ms(imacro,Temperature_ion,dt)
       Enh_gkTlns_m(imacro)= Enh_gkTlns_m(imacro) + gkT * xi_nh_m(imacro)*dt
       Enh_m(imacro)       = Enh_gkTlns_m(imacro) + 0.5d0* Qnh* xi_nh_m(imacro)**2
    endif

    velocity_m(:,:,imacro) = velocity(:,:)
    Rion_m(:,:,imacro)     = Rion(:,:)

    call take_back_mol_into_ucell_ms(imacro)
   
  end subroutine

  subroutine dt_evolve_MD_2_MS(it,imacro)
    use md_ground_state
    implicit none
    integer :: it,imacro, ia
    real(8) :: dt_h,Ework, Temperature_ion, mass_au

    dt_h = dt/2d0
    velocity(:,:) = velocity_m(:,:,imacro)
    force(:,:)    = force_m(:,:,imacro)

    !update ion velocity with dt/2
    Ework = 0d0
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
      !velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h  &
       &                               - fric * velocity(:,ia) * dt_h
    enddo

    !NHC act on velocity with dt/2
    if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
       call apply_nose_hoover_velocity_ms(imacro,dt_h)
    endif

    if (stop_system_momt=='y') call remove_system_momentum(0)
    call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)

    data_local_MD(1,imacro,it) = Temperature_ion
    data_local_MD(2,imacro,it) = Eall_ms(imacro)  !=Uene

    Tene_ms(imacro) = Tion
    Eall_ms(imacro) = Eall_ms(imacro) + Tion

    data_local_MD(3,imacro,it) = Eall_ms(imacro)  !=Uene+Tene

    if(ensemble=="NVT".and.thermostat=="nose-hoover")then
       Enh_m_t(it,imacro)  = Enh_m(imacro)
       Hnvt_m_t(it,imacro) = Enh_m(imacro) + Eall_ms(imacro)
    endif

    velocity_m(:,:,imacro) = velocity(:,:)

  end subroutine

  subroutine init_some_var_ms(imacro)
    use md_ground_state, only: cal_Tion_Temperature_ion
    use force_field_CRK, only: Qai
    implicit none
    integer :: ia,imacro,imol
    real(8) :: Tion, Temperature_ion

    do imol=1,nmol
       Qai_old_ms(:,imol,imacro) = Qai(:,imol)
    enddo

    !!hoge when restart option, this must not be done..?
    velocity(:,:) = velocity_m(:,:,imacro)
    call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
    Eall0_m(imacro) = Uene + Tion

  end subroutine

  subroutine cal_force_energy_CRK_MS(imacro)
    use force_field_CRK, only: Qai,force_intra,Uene_intra,Rxyz_wX,cal_force_energy_CRK,natom_wX_mol
    implicit none
    integer :: ixyz, ia,imacro,imol,ia_wX, ipt, k, iflag_update_Rxyz
    real(8) :: Eti(10) !, Qai_intpl, Qai_hist(10)

    !electric field used in CRK force field
    if(iter2==1) then
       do ixyz=1,3
          k=0
          do ipt = npt_intpl_E-1,1,-1
             k=k+1
             Eti(k) = E_hist_m(ixyz,imacro,ipt)
          enddo
          k=k+1
          Eti(k) = E_m(ixyz,imacro)
          call interpolate(Et(ixyz),Eti,npt_intpl_E)
       enddo
    else
       Et(:)= E_m(:,imacro)
    endif
   !if(imacro==1)write(*,*)"aho",iter2,real(Et(2)),real(Et(3))  !for check
    E_m_iter(:,imacro) = Et(:)


    if(iter2==1) then
       iflag_update_Rxyz = 0
    else 
       iflag_update_Rxyz = 1  !flag to skip intramolecular interaction
       force_intra(:,:) = force_intra_m(:,:,imacro)
       Uene_intra = Uene_intra_ms(imacro)
    endif
    Rion(:,:) = Rion_m(:,:,imacro)

    call cal_force_energy_CRK(iflag_update_Rxyz)

    if(iter2==1) then
       force_intra_m(:,:,imacro) = force_intra(:,:)
       Uene_intra_ms(imacro) = Uene_intra
    endif

    do ia=1,NI
       force_m(:,ia,imacro) = force(:,ia)
    enddo

    do imol=1,nmol
       Qai_old_ms(:,imol,imacro)  = Qai_ms(:,imol,imacro)
    enddo

    do imol=1,nmol
       Qai_ms(:,imol,imacro) = Qai(:,imol)
       Rxyz_wX_ms(:,:,imol,imacro)= Rxyz_wX(:,:,imol)
    enddo
    Uene_ms(imacro) = Uene
    Eall_ms(imacro) = Uene  !electric field energy is included already(need check)

  end subroutine

  subroutine cal_dipole_moment_ms(dipole_mom,dt,imacro)
    ! dipole_mom : (in & out)
    use force_field_CRK, only: Qai, Rxyz_wX, natom_wX_mol, Vuc
    implicit none
    integer :: imol,ia_wX,imacro
    real(8) :: dipole_mom(3), dt
    
    !dipole moment in unit cell
    dipole_mom(:) = 0d0
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       dipole_mom(:) = dipole_mom(:) &
                     + Qai_ms(ia_wX,imol,imacro)*Rxyz_wX_ms(:,ia_wX,imol,imacro)
    enddo
    enddo
    dipole_mom(:) = dipole_mom(:)/aLxyz

  end subroutine

  subroutine cal_dipole_moment_current_ms_EB(dipole_mom,current_matter,dt,imacro)
    ! dipole_mom : (in & out)
    ! claculate Jm(t+dt/2) for use of E(t),B(t)
    use force_field_CRK, only: Qai,Rxyz_wX,natom_wX_mol,Vuc
    implicit none
    integer :: imol,ia_wX,imacro, ii,ia
    real(8) :: dipole_mom(3), current_matter(3), dt
    real(8) :: velocity_wX(3,4,nmol), dQaidt(4,nmol), Qai_mid(4,nmol)

    do imol=1,nmol
       dQaidt(:,imol) = (Qai_ms(:,imol,imacro) - Qai_bk_ms(:,imol,imacro))/dt
       Qai_mid(:,imol)= (Qai_ms(:,imol,imacro) + Qai_bk_ms(:,imol,imacro))*0.5d0
    enddo


    velocity(:,:) = velocity_m_mid(:,:,imacro)
    Rion(:,:)     = Rion_m_mid(:,:,imacro)
    call cal_vel_siteX(velocity_wX)

    !dipole moment in unit cell
    dipole_mom(:)    =0d0
    current_matter(:)=0d0
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       dipole_mom(:) = dipole_mom(:) &
                     + Qai_mid(ia_wX,imol)*Rxyz_wX_ms_mid(:,ia_wX,imol,imacro)
       current_matter(:) = current_matter(:) &
                     -dQaidt(ia_wX,imol) * Rxyz_wX_ms_mid(:,ia_wX,imol,imacro) &
                     - Qai_mid(ia_wX,imol) * velocity_wX(:,ia_wX,imol)
    enddo
    enddo
    dipole_mom(:)     = dipole_mom(:)/aLxyz
    current_matter(:) = current_matter(:)/aLxyz

    !matter current defined by positive charge-->minus sign (less accurate/rough)
    !current_matter(:) = -(dipole_mom(:) - dipole_mom_last(:))/dt

    if(flag_method_newton_Et) then
       do imol=1,nmol
          ii = sum(natom_wX_mol(1:imol)) - natom_wX_mol(imol)
          do ia_wX=1,natom_wX_mol(imol)
             ia = ii+ia_wX
             vpr(:,ia,imacro) = 0.5d0*velocity_wX(:,ia_wX,imol) + Rxyz_wX_ms_mid(:,ia_wX,imol,imacro)/dt
             vmr(:,ia,imacro) = 0.5d0*velocity_wX(:,ia_wX,imol) - Rxyz_wX_ms_mid(:,ia_wX,imol,imacro)/dt
             rvec(:,ia,imacro)  = Rxyz_wX_ms(:,ia_wX,imol,imacro)
             Qvec_bk(ia,imacro) = Qai_bk_ms(ia_wX,imol,imacro)
          enddo
       enddo
    endif

  end subroutine

  subroutine  cal_vel_siteX(velocity_wX)
    use force_field_CRK, only: dist_OX_wat
    implicit none
    real(8) :: velocity_wX(3,4,nmol), vel_O(3)
    integer :: i,j, iO,iH1,iH2,iX1,iX2, iatom_O,iatom_H1,iatom_H2, imol
    real(8) :: drOH1_x_drOH2(3), drOH1(3), drOH2(3), rOH1, rOH2
    real(8) :: e_OH1_x_OH2(3), rOH1rOH2, rOH1rOH2_inv, r2_OH1,r2_OH2
    real(8) :: drOH1_x_drOH2_inv, drOH1_x_drOH2_inv2
    real(8) :: dot_drOH1(3), dot_drOH2(3)
    real(8) :: dot_e_OH1_x_OH2(3), dot_drOH1_x_drOH2(3), dot_drOH1_x_drOH2_abs
    
    iO =1 !(O atom)
    iH1=2 !(H1 atom)
    iH2=3 !(H2 atom)

    iX1=1 !(X1 site)
    iX2=2 !(X2 site)

    do imol=1,nmol

       iatom_O  = mol2atom_top(imol)
       iatom_H1 = mol2atom_top(imol) + 1
       iatom_H2 = mol2atom_top(imol) + 2

       drOH1(:) = Rion(:,iatom_H1) - Rion(:,iatom_O)
       drOH2(:) = Rion(:,iatom_H2) - Rion(:,iatom_O)

!       r2_OH1  = sum(drOH1(:)*drOH1(:))
!       r2_OH2  = sum(drOH2(:)*drOH2(:))
!       rOH1    = sqrt(r2_OH1)
!       rOH2    = sqrt(r2_OH2)
!       rOH1rOH2     = rOH1 * rOH2
!       rOH1rOH2_inv = 1d0/rOH1rOH2

       drOH1_x_drOH2(1)  = drOH1(2)*drOH2(3) - drOH1(3)*drOH2(2)
       drOH1_x_drOH2(2)  = drOH1(3)*drOH2(1) - drOH1(1)*drOH2(3)
       drOH1_x_drOH2(3)  = drOH1(1)*drOH2(2) - drOH1(2)*drOH2(1)
       drOH1_x_drOH2_inv = 1d0/sqrt(sum(drOH1_x_drOH2(:)*drOH1_x_drOH2(:)))
       drOH1_x_drOH2_inv2= drOH1_x_drOH2_inv * drOH1_x_drOH2_inv
       e_OH1_x_OH2(:)    = drOH1_x_drOH2(:) * drOH1_x_drOH2_inv

       dot_drOH1(:) = velocity(:,iatom_H1) - velocity(:,iatom_O)
       dot_drOH2(:) = velocity(:,iatom_H2) - velocity(:,iatom_O)
       
       dot_drOH1_x_drOH2(1) = ( drOH1(2)*dot_drOH2(3) + dot_drOH1(2)*drOH2(3) ) &
                            - ( drOH1(3)*dot_drOH2(2) + dot_drOH1(3)*drOH2(2) )
       dot_drOH1_x_drOH2(2) = ( drOH1(3)*dot_drOH2(1) + dot_drOH1(3)*drOH2(1) ) &
                            - ( drOH1(1)*dot_drOH2(3) + dot_drOH1(1)*drOH2(3) )
       dot_drOH1_x_drOH2(3) = ( drOH1(1)*dot_drOH2(2) + dot_drOH1(1)*drOH2(2) ) &
                            - ( drOH1(2)*dot_drOH2(1) + dot_drOH1(2)*drOH2(1) )

       dot_drOH1_x_drOH2_abs = sum(drOH1_x_drOH2(:)*dot_drOH1_x_drOH2(:)) * drOH1_x_drOH2_inv
       dot_e_OH1_x_OH2(:) = dot_drOH1_x_drOH2(:) * drOH1_x_drOH2_inv - drOH1_x_drOH2(:) * dot_drOH1_x_drOH2_abs * drOH1_x_drOH2_inv2

       vel_O(:)              = velocity(:,iatom_O)
       velocity_wX(:,1,imol) = vel_O(:) + dist_OX_wat * dot_e_OH1_x_OH2(:)  ! X1 site
       velocity_wX(:,2,imol) = vel_O(:) - dist_OX_wat * dot_e_OH1_x_OH2(:)  ! X2 site
       velocity_wX(:,3,imol) = velocity(:,iatom_H1)
       velocity_wX(:,4,imol) = velocity(:,iatom_H2)
       
    enddo
    
  end subroutine cal_vel_siteX
  
end subroutine cmd_maxwell_ms

!===========================================================

subroutine save_matrix_for_newton_method(imacro)
  use force_field_CRK, only: NI_wX, Kmat, invJmat
  implicit none
  integer :: ia,imacro
  real(8) :: vec2_tmp(NI_wX), vec3_tmp(NI_wX)

  do ia=1,NI_wX
     vec2_tmp(ia) = sum( Kmat(ia,:)*rvec(2,:,imacro) )
     vec3_tmp(ia) = sum( Kmat(ia,:)*rvec(3,:,imacro) )
  enddo

  do ia=1,NI_wX
     vec2_tmp_m(ia,imacro) = sum( invJmat(ia,:) * vec2_tmp(:) )
     vec3_tmp_m(ia,imacro) = sum( invJmat(ia,:) * vec3_tmp(:) )
  enddo

end subroutine


subroutine interpolate(Xintpl,Xi,npt)
  !Lagrange interpolation
  integer :: npt
  real(8) :: Xintpl, Xi(10),tmp(10) !Et3,Et2,Et1, tmp1,tmp2

  if(     npt==2) then
     tmp(2) = dble((3-1))/dble((2-1)) * Xi(2)
     tmp(1) = dble((3-2))/dble((1-2)) * Xi(1)
     Xintpl = tmp(1) + tmp(2)

  else if(npt==3) then
     tmp(3) = dble( (4-2)*(4-1))/dble( (3-2)*(3-1) ) * Xi(3)
     tmp(2) = dble( (4-3)*(4-1))/dble( (2-3)*(2-1) ) * Xi(2)
     tmp(1) = dble( (4-3)*(4-2))/dble( (1-3)*(1-2) ) * Xi(1)
     Xintpl = tmp(1) + tmp(2) + tmp(3)

  else if(npt==4) then
     tmp(4) = dble( (5-3)*(5-2)*(5-1))/dble( (4-3)*(4-2)*(4-1) ) * Xi(4)
     tmp(3) = dble( (5-4)*(5-2)*(5-1))/dble( (3-4)*(3-2)*(3-1) ) * Xi(3)
     tmp(2) = dble( (5-4)*(5-3)*(5-1))/dble( (2-4)*(2-3)*(2-1) ) * Xi(2)
     tmp(1) = dble( (5-4)*(5-3)*(5-2))/dble( (1-4)*(1-3)*(1-2) ) * Xi(1)
     Xintpl = tmp(1) + tmp(2) + tmp(3) + tmp(4)

  endif

end subroutine interpolate

!===========================================================
subroutine dt_evolve_EB_1d_cmd
  use Global_variables
  use salmon_communication
  implicit none
  integer :: ix_m,iy_m,iz_m,ix,iy,iz,imacro
  real(8) :: coef1E,coef2E,coef1B,coef2B
  real(8) :: g(2:3), dgdE(2:3,2:3), inv_dgdE(2:3,2:3), dEn(2:3)
  real(8) :: vmr2Q, vmr3Q, det, E_m_new_tmp(3,nmacro)

  iy_m = ny_origin_m
  iz_m = nz_origin_m

  if(flag_use_light_source) call apply_light_source_1d_J

  coef1E = c_light*dt /HX_m
  coef2E = 4*pi*dt
  coef1B = c_light*dt /HX_m

!$omp parallel do default(shared) private(ix_m)
  do ix_m = nx1_m, nx2_m
     E_ms(2,ix_m,iy_m,iz_m) =  E_ms(2,ix_m,iy_m,iz_m) - coef1E * (B_ms(3,ix_m,iy_m,iz_m)-B_ms(3,ix_m-1,iy_m,iz_m)) &
                          & + Jm_ms(2,ix_m,iy_m,iz_m) * coef2E
     E_ms(3,ix_m,iy_m,iz_m) =  E_ms(3,ix_m,iy_m,iz_m) + coef1E * (B_ms(2,ix_m,iy_m,iz_m)-B_ms(2,ix_m-1,iy_m,iz_m)) &
                          & + Jm_ms(3,ix_m,iy_m,iz_m) * coef2E
  end do
!$omp end parallel do


  if(flag_method_newton_Et) then
     !Newton method for E(t)

     E_m_new_tmp(:,:) = 0d0
     do imacro = nmacro_s, nmacro_e
        ix = macropoint(1,imacro)
        iy = macropoint(2,imacro)
        iz = macropoint(3,imacro)

        g(2) = E_m_iter(2,imacro) - ( E_last(2,ix,iy,iz)            &
                   - coef1E * (B_ms(3,ix,iy,iz)-B_ms(3,ix-1,iy,iz)) &
                   + Jm_ms(2,ix,iy,iz) * coef2E )
        g(3) = E_m_iter(3,imacro) - ( E_last(3,ix,iy,iz)            &
                   + coef1E * (B_ms(2,ix,iy,iz)-B_ms(2,ix-1,iy,iz)) &
                   + Jm_ms(3,ix,iy,iz) * coef2E )

        vmr2Q = sum( vmr(2,:,imacro) * Qvec_bk(:,imacro) )
        vmr3Q = sum( vmr(3,:,imacro) * Qvec_bk(:,imacro) )

        dgdE(2,2) = sum( vpr(2,:,imacro) * vec2_tmp_m(:,imacro) ) + vmr2Q
        dgdE(2,3) = sum( vpr(2,:,imacro) * vec3_tmp_m(:,imacro) ) + vmr2Q
        dgdE(3,2) = sum( vpr(3,:,imacro) * vec2_tmp_m(:,imacro) ) + vmr3Q
        dgdE(3,3) = sum( vpr(3,:,imacro) * vec3_tmp_m(:,imacro) ) + vmr3Q

        dgdE(2,2) = 1d0 + coef2E * dgdE(2,2)/aLxyz
        dgdE(2,3) =     + coef2E * dgdE(2,3)/aLxyz
        dgdE(3,2) =     + coef2E * dgdE(3,2)/aLxyz
        dgdE(3,3) = 1d0 + coef2E * dgdE(3,3)/aLxyz

        !(inverse matrix)
        det = dgdE(2,2)*dgdE(3,3) - dgdE(2,3)*dgdE(3,2)
        inv_dgdE(2,2) =   dgdE(3,3) / det
        inv_dgdE(3,3) =   dgdE(2,2) / det
        inv_dgdE(2,3) = - dgdE(2,3) / det
        inv_dgdE(3,2) = - dgdE(3,2) / det

        dEn(2) = inv_dgdE(2,2) * g(2) + inv_dgdE(2,3) * g(3)
        dEn(3) = inv_dgdE(3,2) * g(2) + inv_dgdE(3,3) * g(3)

        !(new E field)
        E_m_new_tmp(2,imacro) = E_m_iter(2,imacro) - dEn(2)
        E_m_new_tmp(3,imacro) = E_m_iter(3,imacro) - dEn(3)

     enddo

     call comm_summation(E_m_new_tmp, E_m_iter, 3*nmacro, nproc_group_global)
!$omp parallel do private(imacro,ix,iy,iz)
     do imacro = 1,nmacro
        ix = macropoint(1,imacro)
        iy = macropoint(2,imacro)
        iz = macropoint(3,imacro)
        E_ms(:,ix,iy,iz) = E_m_iter(:,imacro)
     enddo
!$omp end parallel do

  endif  !flag_method_newton_Et


  if(flag_use_absorb_bound) call apply_absorb_bound_1d_E

!$omp parallel do default(shared) private(ix_m)
  do ix_m = nx1_m, nx2_m
     B_ms(2,ix_m,iy_m,iz_m) = B_ms(2,ix_m,iy_m,iz_m) + coef1B*(E_ms(3,ix_m+1,iy_m,iz_m)-E_ms(3,ix_m,iy_m,iz_m))
     B_ms(3,ix_m,iy_m,iz_m) = B_ms(3,ix_m,iy_m,iz_m) - coef1B*(E_ms(2,ix_m+1,iy_m,iz_m)-E_ms(2,ix_m,iy_m,iz_m))
  end do
!$omp end parallel do

  if(flag_use_absorb_bound) call apply_absorb_bound_1d_B

contains 

  subroutine  apply_absorb_bound_1d_E
     if(     ABC_method=="Mur") then
        call  apply_absorb_bound_1d_Mur_E
     else if(ABC_method=="PML") then
        call  apply_absorb_bound_1d_PML_E
     else
        stop
     endif
  end subroutine

  subroutine  apply_absorb_bound_1d_B
     if(     ABC_method=="Mur") then
        call  apply_absorb_bound_1d_Mur_B
     else if(ABC_method=="PML") then
        call  apply_absorb_bound_1d_PML_B
     else
        stop
     endif
  end subroutine

  !----------------------------------------------------------------------------------------
  subroutine  apply_absorb_bound_1d_Mur_E   ! 1-th Mur
    implicit none
    integer :: ix1,ix2, iy, iz
    real(8) :: tmp

    iy = ny_origin_m
    iz = nz_origin_m
    tmp = ( c_light * dt - HX_m ) / ( c_light * dt + HX_m )

    !!right absorbing boundary
    !ix1=nx2_m
    !ix2=nx2_m-1
    !E_ms(:,ix1,iy,iz) = E_last(:,ix2,iy,iz) + tmp* ( E_ms(:,ix2,iy,iz) - E_last(:,ix1,iy,iz))

    !left absorbing boundary
    ix1=nx1_m
    ix2=nx1_m+1
    E_ms(:,ix1,iy,iz) = E_last(:,ix2,iy,iz) + tmp* ( E_ms(:,ix2,iy,iz) - E_last(:,ix1,iy,iz))

  end subroutine apply_absorb_bound_1d_Mur_E

  subroutine  apply_absorb_bound_1d_Mur_B   ! 1-th Mur
    implicit none
    integer :: ix1,ix2, iy, iz
    real(8) :: tmp

    iy = ny_origin_m
    iz = nz_origin_m
    tmp = ( c_light * dt - HX_m ) / ( c_light * dt + HX_m )

    !right absorbing boundary
    ix1=nx2_m
    ix2=nx2_m-1
    B_ms(:,ix1,iy,iz) = B_last(:,ix2,iy,iz) + tmp* ( B_ms(:,ix2,iy,iz) - B_last(:,ix1,iy,iz))

    !!left absorbing boundary
    !ix1=nx1_m
    !ix2=nx1_m+1
    !B_ms(:,ix1,iy,iz) = B_last(:,ix2,iy,iz) + tmp* ( B_ms(:,ix2,iy,iz) - B_last(:,ix1,iy,iz))

  end subroutine apply_absorb_bound_1d_Mur_B

  !----------------------------------------------------------------------------------------
  subroutine  apply_absorb_bound_1d_PML_E
    implicit none
    integer :: ix,iy,iz
    real(8) :: pi2dt, sgmE_pml

    iy = ny_origin_m
    iz = nz_origin_m

    pi2dt = 2d0 * pi * dt

!$omp parallel do default(shared) private(ix,sgmE_pml,coef1E,coef2E)
    do ix = nx1_m, nx1_pml
       sgmE_pml = sgm_max_E_pml * (dble(nx1_pml-ix)/dble(nx_pml))**M_pml
       coef1E = ( 1d0 - pi2dt * sgmE_pml ) / ( 1d0 + pi2dt * sgmE_pml )
       coef2E = c_light * dt / HX_m / ( 1d0 + pi2dt * sgmE_pml )
       E_ms(2,ix,iy,iz) = coef1E * E_last(2,ix,iy,iz) - coef2E * (B_ms(3,ix,iy,iz)-B_ms(3,ix-1,iy,iz))
       E_ms(3,ix,iy,iz) = coef1E * E_last(3,ix,iy,iz) + coef2E * (B_ms(2,ix,iy,iz)-B_ms(2,ix-1,iy,iz))
    end do
!$omp end parallel do

!$omp parallel do default(shared) private(ix,sgmE_pml,coef1E,coef2E)
    do ix = nx2_pml, nx2_m
       sgmE_pml = sgm_max_E_pml * (dble(ix-nx2_pml)/dble(nx_pml))**M_pml
       coef1E = ( 1d0 - pi2dt * sgmE_pml ) / ( 1d0 + pi2dt * sgmE_pml )
       coef2E = c_light * dt / HX_m / ( 1d0 + pi2dt * sgmE_pml )
       E_ms(2,ix,iy,iz) = coef1E * E_last(2,ix,iy,iz) - coef2E * (B_ms(3,ix,iy,iz)-B_ms(3,ix-1,iy,iz))
       E_ms(3,ix,iy,iz) = coef1E * E_last(3,ix,iy,iz) + coef2E * (B_ms(2,ix,iy,iz)-B_ms(2,ix-1,iy,iz))
    end do
!$omp end parallel do

    E_ms(:,nx1_m,iy,iz) = 0d0
   !E_ms(:,nx2_m,iy,iz) = 0d0

  end subroutine

  subroutine  apply_absorb_bound_1d_PML_B
    implicit none
    integer :: ix,iy,iz
    real(8) :: tmp, sgmH_pml

    iy = ny_origin_m
    iz = nz_origin_m

    tmp = c_light**2 * dt / (8d0 * pi)

!$omp parallel do default(shared) private(ix,sgmH_pml,coef1B,coef2B)
    do ix = nx1_m, nx1_pml-1
       sgmH_pml = sgm_max_B_pml * (dble(nx1_pml-(ix+0.5d0))/dble(nx_pml))**M_pml
       coef1B = ( 1d0 - tmp * sgmH_pml ) / ( 1d0 + tmp * sgmH_pml )
       coef2B = c_light * dt / HX_m / ( 1d0 + tmp * sgmH_pml )
       B_ms(2,ix,iy,iz) = coef1B*B_last(2,ix,iy,iz) + coef2B*(E_ms(3,ix+1,iy,iz)-E_ms(3,ix,iy,iz))
       B_ms(3,ix,iy,iz) = coef1B*B_last(3,ix,iy,iz) - coef2B*(E_ms(2,ix+1,iy,iz)-E_ms(2,ix,iy,iz))
    end do
!$omp end parallel do

!$omp parallel do default(shared) private(ix,sgmH_pml,coef1B,coef2B)
    do ix = nx2_pml, nx2_m
       sgmH_pml = sgm_max_B_pml * (dble((ix+0.5d0)-nx2_pml)/dble(nx_pml))**M_pml
       coef1B = ( 1d0 - tmp * sgmH_pml ) / ( 1d0 + tmp * sgmH_pml )
       coef2B = c_light * dt / HX_m / ( 1d0 + tmp * sgmH_pml )
       B_ms(2,ix,iy,iz) = coef1B*B_last(2,ix,iy,iz) + coef2B*(E_ms(3,ix+1,iy,iz)-E_ms(3,ix,iy,iz))
       B_ms(3,ix,iy,iz) = coef1B*B_last(3,ix,iy,iz) - coef2B*(E_ms(2,ix+1,iy,iz)-E_ms(2,ix,iy,iz))
    end do
!$omp end parallel do

   !B_ms(:,nx1_m,iy,iz) = 0d0
    B_ms(:,nx2_m,iy,iz) = 0d0

  end subroutine

end subroutine dt_evolve_EB_1d_cmd

!===========================================================
subroutine apply_light_source_1d_J
  use Global_variables
  implicit none
  integer :: iy,iz
  real(8) :: time, func_inc, E0inc(3), Jm_add(3), t2_delay

  iy = ny_origin_m
  iz = nz_origin_m
  time = (iter_save-1)*dt


  call add_light_source_1d_J(rlaser_int_wcm2_1, epdir_re1, ae_shape1, t1_delay, &
                             pulse_tw1, omega1, time, E0inc, func_inc)
  Jm_add(:) = c_light/(2d0*pi*HX_m) * E0inc(:) * func_inc
  Jm_ms(:,ix_light_source,iy,iz) = Jm_add(:)

  t2_delay = t1_delay + t1_t2
  call add_light_source_1d_J(rlaser_int_wcm2_2, epdir_re2, ae_shape2, t2_delay, &
                             pulse_tw2, omega2, time, E0inc, func_inc)
  Jm_add(:) = c_light/(2d0*pi*HX_m) * E0inc(:) * func_inc
  Jm_ms(:,ix_light_source,iy,iz) = Jm_ms(:,ix_light_source,iy,iz) + Jm_add(:)

end subroutine apply_light_source_1d_J


subroutine  add_light_source_1d_J(rlaser_int_wcm2, epdir_re, ae_shape, t_delay, &
                                  pulse_tw, omega, time, E0inc, func_inc)
  use inputoutput, only: au_time_fs
  implicit none
  real(8) :: rlaser_int_wcm2, epdir_re(3), t_delay, pulse_tw, omega
  real(8) :: E0inc_abs, E0inc(3)
  real(8) :: time,t_s,t_e, env, func_inc, t_edge
  real(8) :: omg_chirp, omg_min,omg_max, omg_region
  character(*) :: ae_shape

  if(rlaser_int_wcm2.le.0d0) then
     E0inc(:)  = 0d0
  else
     E0inc_abs = 5.338d-9*sqrt(rlaser_int_wcm2)   ! electric field in a.u.
     E0inc(:)  = E0inc_abs * epdir_re(:)
  endif

  t_s = t_delay
  t_e = t_delay + pulse_tw
  env = 0d0

  if(ae_shape=="Acos2")then

     if(time.ge.t_s .and. time.le.t_e) env = sin((time-t_s)/pulse_tw*pi)**2
     func_inc = env * cos(omega*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="Acos4")then

     if(time.ge.t_s .and. time.le.t_e) env = sin((time-t_s)/pulse_tw*pi)**4
     func_inc = env * cos(omega*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="CW_edge20fs")then

     t_edge = 20d0 / au_time_fs
     if(time.ge.t_s .and. time.le.t_e) then
        if(time.ge.t_s .and. time.le.t_s + t_edge) then
           env = sin((time-t_s)/(2d0*t_edge)*pi)**2
        else if( time.ge.t_e - t_edge .and. time.le.t_e) then
           env = sin((time-(t_e-2d0*t_edge))/(2d0*t_edge)*pi)**2
        else
           env = 1.0d0
        endif
     endif
     func_inc = env * cos(omega*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="CW_edge100fs")then

     t_edge = 100d0 / au_time_fs
     if(time.ge.t_s .and. time.le.t_e) then
        if(time.ge.t_s .and. time.le.t_s + t_edge) then
           env = sin((time-t_s)/(2d0*t_edge)*pi)**2
        else if( time.ge.t_e - t_edge .and. time.le.t_e) then
           env = sin((time-(t_e-2d0*t_edge))/(2d0*t_edge)*pi)**2
        else
           env = 1.0d0
        endif
     endif
     func_inc = env * cos(omega*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="CW_chirp")then  !edge=100fs/chirp_range=1000cm

     t_edge = 100d0 / au_time_fs
     omg_chirp = 0d0
     if(time.ge.t_s .and. time.le.t_e) then
        omg_region = 1000d0   ![cm-1]
        omg_region = omg_region * 1.239842428713634d-4 / 27.2114d0 ![cm-1]->[Hartree]
        omg_min    = omega - 0.5d0 * omg_region
        omg_chirp  = omg_min + (time - t_s)/pulse_tw * omg_region
        if(time.ge.t_s .and. time.le.t_s + t_edge) then
           env = sin((time-t_s)/(2d0*t_edge)*pi)**2
        else if( time.ge.t_e - t_edge .and. time.le.t_e) then
           env = sin((time-(t_e-2d0*t_edge))/(2d0*t_edge)*pi)**2
        else
           env = 1.0d0
        endif
     endif
     func_inc = env * cos(omg_chirp*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="CW_chirpB")then  !edge=100fs/chirp_range=2000cm

     t_edge = 100d0 / au_time_fs
     omg_chirp = 0d0
     if(time.ge.t_s .and. time.le.t_e) then
        omg_region = 2000d0   ![cm-1]
        omg_region = omg_region * 1.239842428713634d-4 / 27.2114d0 ![cm-1]->[Hartree]
        omg_min    = omega - 0.5d0 * omg_region
        omg_chirp  = omg_min + (time - t_s)/pulse_tw * omg_region
        if(time.ge.t_s .and. time.le.t_s + t_edge) then
           env = sin((time-t_s)/(2d0*t_edge)*pi)**2
        else if( time.ge.t_e - t_edge .and. time.le.t_e) then
           env = sin((time-(t_e-2d0*t_edge))/(2d0*t_edge)*pi)**2
        else
           env = 1.0d0
        endif
     endif
     func_inc = env * cos(omg_chirp*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="Acos2_chirp")then  !chirp_range=1000cm
     omg_region = 1000d0   ![cm-1]
     omg_region = omg_region * 1.239842428713634d-4 / 27.2114d0 ![cm-1]->[Hartree]
     omg_min    = omega - 0.5d0 * omg_region
     omg_chirp  = omg_min + (time - t_s)/pulse_tw * omg_region
     if(time.ge.t_s .and. time.le.t_e) env = sin((time-t_s)/pulse_tw*pi)**2
     func_inc = env * cos(omg_chirp*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="Acos2_chirpB")then  !chirp_range=2000cm
     omg_region = 2000d0   ![cm-1]
     omg_region = omg_region * 1.239842428713634d-4 / 27.2114d0 ![cm-1]->[Hartree]
     omg_min    = omega - 0.5d0 * omg_region
     omg_chirp  = omg_min + (time - t_s)/pulse_tw * omg_region
     if(time.ge.t_s .and. time.le.t_e) env = sin((time-t_s)/pulse_tw*pi)**2
     func_inc = env * cos(omg_chirp*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="Acos2_chirpB2")then  !chirp_range=2000cm (reverse)
     omg_region = 2000d0   ![cm-1]
     omg_region = omg_region * 1.239842428713634d-4 / 27.2114d0 ![cm-1]->[Hartree]
    !omg_min    = omega - 0.5d0 * omg_region
    !omg_chirp  = omg_min + (time - t_s)/pulse_tw * omg_region
     omg_max    = omega + 0.5d0 * omg_region
     omg_chirp  = omg_max - (time - t_s)/pulse_tw * omg_region
     if(time.ge.t_s .and. time.le.t_e) env = sin((time-t_s)/pulse_tw*pi)**2
     func_inc = env * cos(omg_chirp*((time-t_s-pulse_tw*0.5d0)))

  else if(ae_shape=="none")then
     func_inc = 0d0
  else
     func_inc = 0d0
     if(comm_is_root(nproc_id_global)) write(*,*) "#No such name of ae_shape:", trim(ae_shape)
  endif

end subroutine

subroutine  init_absorb_bound_1d_PML
  implicit none
  real(8) :: tmp_pi2dt, tmp_coef

  sgm_max_E_pml = - (log(R_pml)) *(M_pml+1)/(2d0 * HX_m*nx_pml * (4d0*pi/c_light) )
  sgm_max_B_pml = sgm_max_E_pml * (4*pi/c_light)**2

  nx1_pml = nx1_m + nx_pml -1
  nx2_pml = nx2_m - nx_pml +1

  if(comm_is_root(nproc_id_global)) write(*,*) "sgm_max_E_pml= ", real(sgm_max_E_pml)
  tmp_pi2dt = 2d0 * pi * dt
  tmp_coef  = (1d0 - tmp_pi2dt * sgm_max_E_pml)/(1d0 + tmp_pi2dt * sgm_max_E_pml)
  if(tmp_coef.lt.0d0) then
     if(comm_is_root(nproc_id_global)) write(*,*) "Error: coef is negative", real(tmp_coef)
     stop
  endif

end subroutine

subroutine apply_nose_hoover_velocity_ms(imacro,dt_h)
  implicit none
  integer :: ia,imacro
  real(8) :: tmp_exp,dt_h 
  
  tmp_exp = exp(-xi_nh_m(imacro)*dt_h)
  do ia=1,NI
!     velocity_m(:,ia,imacro) = velocity_m(:,ia,imacro) * tmp_exp
     velocity(:,ia) = velocity(:,ia) * tmp_exp
  enddo
  
  return
end subroutine

subroutine apply_nose_hoover_thermostat_ms(imacro,Temp,dt_f)
  implicit none
  integer :: imacro
  real(8) :: Temp,dt_f 

  xi_nh_m(imacro) = xi_nh_m(imacro) + (Temp/temperature0_ion-1d0)/(thermostat_tau**2d0)*dt_f 

  return
end subroutine apply_nose_hoover_thermostat_ms

subroutine calc_energy_joule_EB()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: elec_mid(3), jm_mid(3), ohm_mid
  !$omp parallel do collapse(3) default(shared) private( &
  !$omp     ix_m, iy_m, iz_m, elec_mid, jm_mid, ohm_mid &
  !$omp     )
  do iz_m = nz1_m, nz2_m
  do iy_m = ny1_m, ny2_m
  do ix_m = nx1_m, nx2_m
     !use at t+dt/2
     jm_mid   = Jm_ms(:,ix_m,iy_m,iz_m)
     elec_mid = (E_ms(:,ix_m,iy_m,iz_m) + E_last(:,ix_m,iy_m,iz_m))*0.5d0
     ohm_mid  = sum(-jm_mid * elec_mid)
     energy_joule_ms(ix_m,iy_m,iz_m) = energy_joule_ms(ix_m,iy_m,iz_m) &
                                     & + ohm_mid * aLxyz * dt 
  end do
  end do
  end do
  !$omp end parallel do
  return
end subroutine calc_energy_joule_EB

!===========================================================
subroutine calc_energy_joule_subtract_light_source_EB()
  use Global_variables
  implicit none
  integer :: ix, iy, iz
  real(8) :: elec_mid(3), jm_mid(3), ohm_mid

  iy = ny_origin_m
  iz = nz_origin_m
  ix = ix_light_source
  !use at t+dt/2
  jm_mid   = Jm_ms(:,ix,iy,iz)   
  elec_mid = ( E_ms(:,ix,iy,iz) + E_last(:,ix,iy,iz) )*0.5d0
  ohm_mid  = sum(-jm_mid * elec_mid)
  energy_joule_ms(ix,iy,iz) = energy_joule_ms(ix,iy,iz) &
                              - ohm_mid * aLxyz * dt  !subtract

  return
end subroutine calc_energy_joule_subtract_light_source_EB

!===========================================================
subroutine calc_energy_elemag_EB()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: e2, b2, tmp

  tmp = (1d0 / (8d0 * pi)) * aLxyz
  !$omp parallel do collapse(3) default(shared) private(ix_m,iy_m,iz_m,e2,b2)
  do iz_m = nz1_m, nz2_m
  do iy_m = ny1_m, ny2_m
  do ix_m = nx1_m, nx2_m
    !e2 = sum(E_ms(:,ix_m,iy_m,iz_m)**2)
     e2 = sum(((E_ms(:,ix_m,iy_m,iz_m)+E_last(:,ix_m,iy_m,iz_m))*0.5d0)**2)
     b2 = sum(B_ms(:,ix_m,iy_m,iz_m)**2)
     energy_elemag_ms(ix_m,iy_m,iz_m) = tmp * (e2 + b2)
  end do
  end do
  end do
  !$omp end parallel do
  return
end subroutine calc_energy_elemag_EB
!===========================================================
subroutine calc_total_energy_EB()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: e_em, e_ej, e_ex
  real(8) :: cell_scalling 
  
  e_em=0d0
  e_ej=0d0
  e_ex=0d0
  cell_scalling = HX_m / aLx   !for 1D
  
!$omp parallel do collapse(3) default(shared) private(ix_m,iy_m,iz_m) reduction(+: e_em,e_ej,e_ex)
  do iz_m = nz1_m, nz2_m
  do iy_m = ny1_m, ny2_m
  do ix_m = nx1_m, nx2_m
     e_em = e_em + energy_elemag_ms(ix_m,iy_m,iz_m) * cell_scalling
     e_ej = e_ej + energy_joule_ms(ix_m,iy_m,iz_m)  * cell_scalling
     e_ex = e_ex + energy_elec_ms(ix_m,iy_m,iz_m)   * cell_scalling
  end do
  end do
  end do
!$end omp parallel do
  total_energy_elemag_old = total_energy_elemag
  total_energy_absorb_old = total_energy_absorb
  total_energy_em_old = total_energy_em
  total_energy_elec_old = total_energy_elec
  total_energy_elemag = e_em
  total_energy_absorb = e_ej
  total_energy_em = e_em + e_ej
  total_energy_elec = e_ex
  return
end subroutine calc_total_energy_EB

!===========================================================

subroutine read_ini_coor_vel_ms_cmd
  use salmon_global
  use salmon_communication
  use salmon_parallel
  use global_variables
  use misc_routines
  use inputoutput, only: au_length_aa
  implicit none
  integer :: unit_ini_rv,ia,j, imacro
  character(1024) :: file_ini_rv
  real(8) :: ulength_to_au

  do imacro = nmacro_s, nmacro_e
 
     select case(unit_length)
     case('au','a.u.')
        ulength_to_au   = 1d0
     case('AA','angstrom','Angstrom')
        ulength_to_au   = 1d0/au_length_aa
     end select

     ! file for initial atomic coordinate and velocity in each macro grid
     ! must be "(parent dir)/multiscale/M00XXXX/ini_coor_vel.dat"
     ! the file format is
     !  r_x  r_y  r_z  v_x  v_y  v_z  (for 1th atom)
     !  r_x  r_y  r_z  v_x  v_y  v_z  (for 2th atom)
     !    ..................
     !  r_x  r_y  r_z  v_x  v_y  v_z  (for NI-th atom)
     ! (unit of coordinate is A or a.u., but velocity must be a.u.)
     if(comm_is_root(nproc_id_tdks)) then
        unit_ini_rv = 5000 + imacro
        file_ini_rv = trim(dir_ms_M(imacro))//'ini_coor_vel.dat'
        open(unit_ini_rv,file=trim(file_ini_rv),status="old")
        do ia=1,NI
           read(unit_ini_rv,*) (Rion(j,ia),j=1,3),(velocity(j,ia),j=1,3)
        enddo
        Rion(:,:)    = Rion(:,:)*ulength_to_au
        Rion_eq(:,:) = Rion(:,:)
        close(unit_ini_rv)
     endif
           
     call comm_bcast(Rion,    nproc_group_tdks)
     call comm_bcast(Rion_eq, nproc_group_tdks)
     call comm_bcast(velocity,nproc_group_tdks)

     Rion_m(:,:,imacro)     = Rion(:,:)
     Rion_eq_m(:,:,imacro)  = Rion_m(:,:,imacro)
     velocity_m(:,:,imacro) = velocity(:,:)
           
  enddo !imacro

end subroutine


subroutine read_write_ms_cmd(iflag_read_write_ms)
  use salmon_global
  use salmon_communication
  use salmon_parallel
  use global_variables
  use misc_routines
  use force_field_CRK, only: Qai
  implicit none
  integer,intent(in) :: iflag_read_write_ms
  integer :: imacro
  integer :: nfile_md_ms, nfile_ae_ms, nfile_other_ms
  integer :: ix_m,iy_m,iz_m
  character(1024) :: md_file_ms, ae_file_ms, dir_ae_file_ms, other_file_ms
  real(8),allocatable :: energy_joule_ms_tmp(:,:,:)

  !copy from GS/io_gs_wfn_k.f90
  character(256) :: rt_wfn_directory
  integer,parameter :: nfile_md     = 43
  integer,parameter :: nfile_ae     = 44
  integer,parameter :: nfile_al     = 45
  integer,parameter :: iflag_read_rt = 0
  integer,parameter :: iflag_write_rt= 1
  !--------------------------------------

  if(.not. allocated(energy_joule_ms_tmp)) &
     allocate( energy_joule_ms_tmp(nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m) )
  if(.not. allocated(D1_pml)) then
     allocate( D1_pml(1:3,mx1_m:mx2_m,my1_m:my2_m,mz1_m:mz2_m) )
     allocate( D2_pml(1:3,mx1_m:mx2_m,my1_m:my2_m,mz1_m:mz2_m) )
  endif

  if (comm_is_root(nproc_id_global)) then
     nfile_ae_ms = 7000
     dir_ae_file_ms = trim(dir_ms)//'rt_ae_field/'
     if(iflag_read_write_ms==iflag_write_rt) call create_directory(dir_ae_file_ms)
     ae_file_ms = trim(dir_ae_file_ms)//'ae_field'
     open(nfile_ae_ms,file=trim(ae_file_ms),form='unformatted')
     select case(iflag_read_write_ms)
     case(iflag_write_rt)
        write(nfile_ae_ms) Ac_new_ms,Ac_ms,Ac_old_ms, &
                           Jm_new_ms,Jm_ms,D1_pml, D2_pml
     case(iflag_read_rt )
        read(nfile_ae_ms ) add_Ac_new_ms,add_Ac_ms,Ac_old_ms, &
                           Jm_new_ms,Jm_ms, D1_pml, D2_pml
     end select
     close(nfile_ae_ms)
  end if

  do imacro = nmacro_s, nmacro_e

     write (rt_wfn_directory,'(A,A)') trim(dir_ms_M(imacro)),'/rt_wfn_k/'

     nfile_md_ms     = 7000 + nproc_id_global
     nfile_other_ms  = nfile_md_ms

     if(comm_is_root(nproc_id_tdks)) then
        if(iflag_read_write_ms==iflag_write_rt) call create_directory(rt_wfn_directory)
        md_file_ms = trim(rt_wfn_directory)//'Rion_velocity'
        open(nfile_md_ms,file=trim(md_file_ms),form='unformatted')
        select case(iflag_read_write_ms)
        case(iflag_write_rt); write(nfile_md_ms) Rion_m(:,:,imacro),velocity_m(:,:,imacro), Qai_ms(:,:,imacro)
        case(iflag_read_rt ); read(nfile_md_ms ) Rion,velocity,Qai
        end select
        close(nfile_md)

        other_file_ms = trim(rt_wfn_directory)//'others'
        open(nfile_other_ms,file=trim(other_file_ms),form='unformatted')
        ix_m = macropoint(1,imacro)
        iy_m = macropoint(2,imacro)
        iz_m = macropoint(3,imacro)
        select case(iflag_read_write_ms)
        case(iflag_write_rt); write(nfile_other_ms) Eall0_m(imacro), energy_joule_ms(ix_m,iy_m,iz_m)
        case(iflag_read_rt ); read(nfile_other_ms ) Eall0_m(imacro), energy_joule_ms_tmp(ix_m,iy_m,iz_m)
        end select
        close(nfile_other_ms)

     end if

     call comm_sync_all


     if(iflag_read_write_ms == iflag_read_rt) then

        call comm_bcast(Rion,    nproc_group_tdks)
        call comm_bcast(velocity,nproc_group_tdks)
        call comm_bcast(Qai,     nproc_group_tdks)
        Rion_eq(:,:)           = Rion(:,:)
        Rion_m(:,:,imacro)     = Rion(:,:)
        Rion_eq_m(:,:,imacro)  = Rion_m(:,:,imacro)
        velocity_m(:,:,imacro) = velocity(:,:)
        Qai_ms(:,:,imacro)     = Qai(:,:)

        call comm_summation(energy_joule_ms_tmp,energy_joule_ms,&
          & (nx2_m-nx1_m+1)*(ny2_m-ny1_m+1)*(nz2_m-nz1_m+1),nproc_group_global)
        call comm_bcast(add_Ac_ms,     nproc_group_global)
        call comm_bcast(add_Ac_new_ms, nproc_group_global)
        call comm_bcast(Ac_old_ms,     nproc_group_global)
        call comm_bcast(Jm_new_ms,     nproc_group_global)
        call comm_bcast(Jm_ms,         nproc_group_global)
        call comm_bcast(D1_pml,        nproc_group_global)
        call comm_bcast(D2_pml,        nproc_group_global)
        
     endif

  enddo  !imacro

  end subroutine read_write_ms_cmd

end module

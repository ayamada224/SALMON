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
!This file is "force_field_CRK.f90"
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module force_field_CRK
  use Global_Variables
  use salmon_parallel
  use salmon_communication
  implicit none
  logical :: flag_use_over_cell, flag_use_Et
  integer :: NI_wX, iter_max,ncell_max(3),nc_max
  integer,allocatable :: natom_wX_mol_s(:), natom_wX_mol(:)
  real(8) :: dQave_thresh
  real(8),allocatable :: Rxyz_wX(:,:,:)
  !each molecule
  integer,allocatable :: iatom_center(:)
  integer,allocatable :: nbnd_r(:), iatom1_bnd_r(:,:), iatom2_bnd_r(:,:)
  integer,allocatable :: nbnd_a(:), iatom1_bnd_a(:,:), iatom2_bnd_a(:,:)
  real(8),allocatable :: r0_bnd_r(:,:), r0_bnd_a(:,:), r0_bnd_c(:,:)
  real(8),allocatable :: kbnd_r(:,:,:,:), kbnd_a(:,:)
  real(8),allocatable :: kbnd_ra(:,:,:), kbnd_rr(:,:,:), kbnd_c(:,:)
  real(8),allocatable :: charge_cln(:,:), epsilon_vdw(:,:), sigma_vdw(:,:)
  real(8) :: the0_HOH, dist_OX_wat, Qeq_wat(4)
  real(8) :: dQdS_wat(4,3), Keq_wat(4,4), dKdS_wat(4,4,3)
  !system
  integer :: nb_r, nb_a, nb_rr, nb_ra
  integer,allocatable :: iatom1_b_r(:), iatom2_b_r(:)
  integer,allocatable :: iatom1_b_a(:), iatom2_b_a(:)
  integer,allocatable :: ibond1_b_rr(:), ibond2_b_rr(:)
  integer,allocatable :: ibond_b_ra(:),  iangle_b_ra(:)
  real(8),allocatable :: r0_b_r(:), r0_b_a(:), r0_b_c(:)
  real(8),allocatable :: kb_r(:,:,:), kb_a(:), kb_ra(:), kb_rr(:), kb_c(:)
  real(8),allocatable :: Qeq(:), Qai(:,:), eps_vdw(:), sgm_vdw(:)
  !switching function for dump_f
  real(8) :: R1_switch,R2_switch, R1_switch2,R2_switch2, R2_R1_swirtch3_inv
  real(8) :: xi_gauss_cln, gm_gauss_cln, coef_derf_gm
  !ewald
  integer :: nG_cmd
  real(8) :: r_cutoff,G_cutoff,Vuc, a_ewald,coef_derfc,coef_derf1
  real(8),allocatable :: Gvec(:,:),G2(:), sumQicosGri(:),sumQisinGri(:)

contains

  subroutine allocate_ff_CRK
    implicit none
    integer :: max_nbnd_r, max_nbnd_a, max_natom

    !change the maximum value if you add larger molecule
    max_natom  = 5
    max_nbnd_r = 5
    max_nbnd_a = 5

    allocate( natom_mol_s(nmol_s) )    ! # of atoms of the molecule species
    allocate( natom_wX_mol_s(nmol_s) ) ! # of atoms of the molecule species including X sites

    allocate( iatom_center(nmol_s) )
    allocate( nbnd_r(nmol_s), nbnd_a(nmol_s) )
    allocate( iatom1_bnd_r(nmol_s,max_nbnd_r), iatom2_bnd_r(nmol_s,max_nbnd_r) )
    allocate( iatom1_bnd_a(nmol_s,max_nbnd_a), iatom2_bnd_a(nmol_s,max_nbnd_a) )
    allocate( r0_bnd_r(nmol_s,max_nbnd_r), r0_bnd_a(nmol_s,max_nbnd_a) )
    allocate( r0_bnd_c(nmol_s,max_nbnd_r) )
    allocate( kbnd_r(nmol_s,max_nbnd_r,2:6,2) )
    allocate( kbnd_a(nmol_s,max_nbnd_a) )
    allocate( kbnd_ra(nmol_s,max_nbnd_r,max_nbnd_a) )
    allocate( kbnd_rr(nmol_s,max_nbnd_r,max_nbnd_r) )
    allocate( kbnd_c(nmol_s,max_nbnd_r) )

    allocate( charge_cln(nmol_s,max_natom) )
    allocate( epsilon_vdw(nmol_s,max_natom), sigma_vdw(nmol_s,max_natom) )

    iatom1_bnd_r(:,:)=0 ; iatom2_bnd_r(:,:)=0
    iatom1_bnd_a(:,:)=0 ; iatom2_bnd_a(:,:)=0
    r0_bnd_r(:,:)=0d0   ; r0_bnd_a(:,:)=0d0    ; r0_bnd_c(:,:)=0d0
    kbnd_r(:,:,:,:)=0d0 ; kbnd_a(:,:)=0d0
    kbnd_ra(:,:,:)=0d0  ; kbnd_rr(:,:,:)=0d0   ; kbnd_c(:,:)=0d0
    charge_cln(:,:)=0d0 ; epsilon_vdw(:,:)=0d0 ; sigma_vdw(:,:)=0d0

  end subroutine

  subroutine load_force_field_CRK(is,name_mol)
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer  is, ib_r,jb_r,ib_a
    real(8) :: Navogadro, eV2J, au2kJmol
    parameter( Navogadro = 6.022140857d23 )
    parameter( eV2J = 1.60218d-19 )
    character(*) :: name_mol

    if(name_mol=="WAT"   .or. &
       name_mol=="WATER" .or. &
       name_mol=="Water" .or. &
       name_mol=="water" ) then
       ! atom order: O, H1, H2

       natom_mol_s(is)  = 3  !O, H1, H2
       iatom_center(is) = 1  !Oxigen

       natom_wX_mol_s(is) = 4 !X1, X2, H1, H2 (only for coulomb interaction)

       nbnd_r(is) = 2
       nbnd_a(is) = 1

       ! r12 for O-H1 and O-H2 (1-2 interaction)
       do ib_r=1,nbnd_r(is)
          iatom1_bnd_r(is,ib_r) = 1       !O
          iatom2_bnd_r(is,ib_r) = ib_r+1  !H1 or H2

          r0_bnd_r(is,ib_r) = 0.9572d0/au_length_aa ![A]->[au]

          kbnd_r(is,ib_r,2,1:2) = (/  0.269d0,  0.269d0 /)  !n=2 [au]
          kbnd_r(is,ib_r,3,1:2) = (/ -0.450d0, -0.500d0 /)  !n=3 [au]
          kbnd_r(is,ib_r,4,1:2) = (/  0.550d0,  3.000d0 /)  !n=4 [au]
          kbnd_r(is,ib_r,5,1:2) = (/ -3.000d0, -3.900d0 /)  !n=5 [au]
          kbnd_r(is,ib_r,6,1:2) = (/  7.000d0, 20.000d0 /)  !n=6 [au]

          !for prevent polarization catastrophy (see ishiyama's paper in 2011)
          r0_bnd_c(is,ib_r) = 0.16d0/au_length_aa ![A]->[au]
          kbnd_c(is,ib_r)   = 5d2  ! [au]
          !(following is suggestion from ishiyama-kun)
          !r0_bnd_c(is,ib_r) = 0.15d0/au_length_aa ![A]->[au]
          !kbnd_c(is,ib_r)   = 2d3  ! [au]
       enddo

       ! r13 for H-H (1-3 interaction)
       ib_a=1
       iatom1_bnd_a(is,ib_a) = 2  !H1
       iatom2_bnd_a(is,ib_a) = 3  !H2
       the0_HOH = 104.52d0 * pi/180d0 ! H-O-H = 104.52 [deg]
       r0_bnd_a(is,ib_a) = 2d0*(0.9572d0/au_length_aa)*sin(the0_HOH/2d0)  ![au]
       kbnd_a(is,ib_a)   = 0.109d0 ![au]
       !write(*,*) "r0(H-H)=", r0_bnd_a(is,ib_a)

       ! r12*r13
       ib_a=1
       do ib_r=1,nbnd_r(is)
          kbnd_ra(is,ib_r,ib_a)  = -0.0822  ![au]
       enddo

      ! r12*r12
       do ib_r=1,nbnd_r(is)
       do jb_r=ib_r+1,nbnd_r(is)
          kbnd_rr(is,ib_r,jb_r) = 0.0678d0  ![au]
       enddo
       enddo

       dist_OX_wat = 0.55d0 / au_length_aa ![bohr] for X site

       !these charges are not used actually !just for test
       charge_cln(is,1) =-0.329d0*2  ![e]
       charge_cln(is,2) = 0.329d0
       charge_cln(is,3) = 0.329d0

       !!just test
       !charge_cln(is,1) =-0.4d0*2  ![e]
       !charge_cln(is,2) = 0.4d0
       !charge_cln(is,3) = 0.4d0

       !these charges are used
       Qeq_wat(1) = -0.329d0  ![e] X1 site
       Qeq_wat(2) = -0.329d0  ![e] X2 site
       Qeq_wat(3) =  0.329d0  ![e] H1 atoms
       Qeq_wat(4) =  0.329d0  ![e] H2 atoms

       !parameter for conformation-dependent partial charge [a.u.]
       !                 (/   X1, X2, H1, H2  /)
       dQdS_wat(1:4,1) = (/  0.097d0,  0.097d0, -0.097d0, -0.097d0 /) !for S1(t=1)
       dQdS_wat(1:4,2) = (/ -0.060d0, -0.060d0,  0.060d0,  0.060d0 /) !for S2(t=2)
       dQdS_wat(1:4,3) = (/  0.0d0,    0.0d0,   -0.049d0,  0.049d0 /) !for S3(t=3)

       !parameter of Keq_abi
       !                 (/   X1, X2, H1, H2  /)
       Keq_wat(1,1:4)  = (/ -4.298d0,  0.248d0,  2.025d0,  2.025d0 /)  !a=X1
       Keq_wat(2,1:4)  = (/  0.248d0, -4.298d0,  2.025d0,  2.025d0 /)  !a=X2
       Keq_wat(3,1:4)  = (/  2.025d0,  2.025d0, -3.271d0, -0.779d0 /)  !a=H1
       Keq_wat(4,1:4)  = (/  2.025d0,  2.025d0, -0.779d0, -3.271d0 /)  !a=H2

       !parameter of dKab/dSt
       !                 (/   X1, X2, H1, H2  /)
       !(for t=1 : S1)
       dKdS_wat(1,1:4,1)  = (/  0.034d0,  1.120d0, -0.579d0, -0.579d0 /)  !a=X1
       dKdS_wat(2,1:4,1)  = (/  1.120d0,  0.034d0, -0.579d0, -0.579d0 /)  !a=X2
       dKdS_wat(3,1:4,1)  = (/ -0.579d0, -0.579d0,  0.621d0,  0.537d0 /)  !a=H1
       dKdS_wat(4,1:4,1)  = (/ -0.579d0, -0.579d0,  0.537d0,  0.621d0 /)  !a=H2
       !(for t=2 : S2)
       dKdS_wat(1,1:4,2)  = (/ -2.706d0, -2.440d0,  2.573d0,  2.573d0 /)  !a=X1
       dKdS_wat(2,1:4,2)  = (/ -2.440d0, -2.706d0,  2.573d0,  2.573d0 /)  !a=X2
       dKdS_wat(3,1:4,2)  = (/  2.573d0,  2.573d0, -1.784d0, -3.362d0 /)  !a=H1
       dKdS_wat(4,1:4,2)  = (/  2.573d0,  2.573d0, -3.362d0, -1.784d0 /)  !a=H1
       !(for t=3 : S3)
       dKdS_wat(1,1:4,3)  = (/  0.000d0,  0.000d0, -0.671d0,  0.671d0 /)  !a=X1
       dKdS_wat(2,1:4,3)  = (/  0.000d0,  0.000d0, -0.671d0,  0.671d0 /)  !a=X2
       dKdS_wat(3,1:4,3)  = (/ -0.671d0, -0.671d0,  1.342d0,  0.000d0 /)  !a=H1
       dKdS_wat(4,1:4,3)  = (/  0.671d0,  0.671d0,  0.000d0, -1.342d0 /)  !a=H2

       !VDW
       au2kJmol = au_energy_ev * eV2J / 1d3 * Navogadro

       epsilon_vdw(is,1) = 0.6496d0/au2kJmol  ![kJ/mol]->[au]
       epsilon_vdw(is,2) = 0d0
       epsilon_vdw(is,3) = 0d0

       sigma_vdw(is,1) = 3.205d0/au_length_aa  ![A]->[au]
       sigma_vdw(is,2) = 0d0
       sigma_vdw(is,3) = 0d0


    !else if(name_mol=="xxxxx") then
    !   ! Dammy molecule
    !   natom_mol_s(is) = 1
    
    else
       call Err_finalize("no such molecular name")
       stop
    endif

  end subroutine

  subroutine prepare_force_field_CRK
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer :: imol, is, iatom, ia, is_wat
    integer :: ib_r, jb_r, ib_a, iib_r, iib_a, iib_rr, iib_ra, ib_r_last_mol, ib_a_last_mol
    real(8) :: Q0ai(4,nmol), dQdr_wat(3,4,3,nmol)
    real(8) :: Kai(4,4,nmol), dKdr_wat(3,4,4,3,nmol) 

    flag_use_Et = .true.

    !(calculate nb_r, nb_a)
    nb_r=0
    nb_a=0
    nb_rr=0
    nb_ra=0
    do imol=1,nmol
       is = mol2species(imol)
       mol2atom_cnt(imol) = mol2atom_top(imol)-1 + iatom_center(is)  !center atom
       do ib_r=1,nbnd_r(is)
          nb_r = nb_r + 1
          do jb_r=ib_r+1,nbnd_r(is)
             nb_rr = nb_rr + 1
          enddo
          do ib_a=1,nbnd_a(is)
             nb_ra = nb_ra + 1
          enddo
       enddo
       do ib_a=1,nbnd_a(is)
          nb_a = nb_a + 1
       enddo
    enddo

    !(allocate)
    allocate( Rxyz_wX(3,4,nmol) )
    allocate( natom_wX_mol(nmol) )
    allocate( iatom1_b_r(nb_r), iatom2_b_r(nb_r) )
    allocate( iatom1_b_a(nb_a), iatom2_b_a(nb_a) )
    allocate( r0_b_r(nb_r), r0_b_a(nb_a), r0_b_c(nb_r) )
    allocate( kb_r(nb_r,2:6,2), kb_a(nb_a) )
    allocate( kb_ra(nb_ra), kb_rr(nb_rr), kb_c(nb_r) )
    allocate( ibond1_b_rr(nb_rr), ibond2_b_rr(nb_rr) )
    allocate( ibond_b_ra(nb_ra),  iangle_b_ra(nb_ra) )
    allocate( Qeq(NI), eps_vdw(NI), sgm_vdw(NI) )
    allocate( Qai(4,nmol) )

    !(set parameters)
    iib_r  = 0
    iib_a  = 0
    iib_rr = 0
    iib_ra = 0
    NI_wX  = 0
    ib_r_last_mol =0
    ib_a_last_mol =0
    do imol=1,nmol
       is = mol2species(imol)
       natom_wX_mol(imol) = natom_wX_mol_s(is)
       NI_wX = NI_wX + natom_wX_mol(imol)
       do ib_r=1,nbnd_r(is)
          iib_r = iib_r + 1
          iatom1_b_r(iib_r) = iatom1_bnd_r(is,ib_r) + mol2atom_top(imol)-1
          iatom2_b_r(iib_r) = iatom2_bnd_r(is,ib_r) + mol2atom_top(imol)-1
          r0_b_r(iib_r) = r0_bnd_r(is,ib_r)
          kb_r(iib_r,2:6,:) = kbnd_r(is,ib_r,2:6,:)

          r0_b_c(iib_r) = r0_bnd_c(is,ib_r)
          kb_c(iib_r)   = kbnd_c(is,ib_r)
          do jb_r=ib_r+1,nbnd_r(is)
             iib_rr = iib_rr + 1
             kb_rr(iib_rr) = kbnd_rr(is,ib_r,jb_r)
             ibond1_b_rr(iib_rr) = ib_r + ib_r_last_mol
             ibond2_b_rr(iib_rr) = jb_r + ib_r_last_mol
          enddo
          do ib_a=1,nbnd_a(is)
             iib_ra = iib_ra + 1
             kb_ra(iib_ra) = kbnd_ra(is,ib_r,ib_a)
             ibond_b_ra(iib_ra)  = ib_r + ib_r_last_mol
             iangle_b_ra(iib_ra) = ib_a + ib_a_last_mol
          enddo
       enddo

       do ib_a=1,nbnd_a(is)
          iib_a = iib_a + 1
          iatom1_b_a(iib_a) = iatom1_bnd_a(is,ib_a) + mol2atom_top(imol)-1
          iatom2_b_a(iib_a) = iatom2_bnd_a(is,ib_a) + mol2atom_top(imol)-1
          r0_b_a(iib_a) = r0_bnd_a(is,ib_a)
          kb_a(iib_a) = kbnd_a(is,ib_a)
       enddo

       do iatom=1,natom_mol(imol)
          ia = mol2atom_top(imol) + iatom-1
          Qeq(ia)     = charge_cln(is,iatom)
          eps_vdw(ia) = epsilon_vdw(is,iatom)
          sgm_vdw(ia) = sigma_vdw(is,iatom)
       enddo

       ib_r_last_mol = ib_r_last_mol + nbnd_r(is)
       ib_a_last_mol = ib_a_last_mol + nbnd_a(is)

    enddo  ! imol

    !(some parameters): read from input file (cmd_sc.inp)
    !iter_max = 50
    !dQave_thresh = 1d-7 ![au]
    !r_cutoff = 13.0d0/au_length_aa ! [Bohr]
    !G_cutoff = 1.47d0*au_length_aa ! [1/Bohr]

    !(cutoff for real space)
    flag_use_over_cell=.false.
    if(r_cutoff .gt. minval(aL)/2d0) then
      !r_cutoff = minval(aL)/2d0
       flag_use_over_cell=.true.
       ncell_max(:) = int(r_cutoff/aL(:))+1
    endif
    if(flag_use_over_cell) then
       nc_max = (2*ncell_max(1)+1)*(2*ncell_max(2)+1)*(2*ncell_max(3)+1)
    else
       nc_max = 1
    endif
    if(comm_is_root(nproc_id_global)) write(*,*) "flag_use_over_cell=",flag_use_over_cell
    if(comm_is_root(nproc_id_global)) write(*,*) "R-cutoff radius =",r_cutoff*au_length_aa," [A]"
    if(comm_is_root(nproc_id_global)) write(*,*) "G-cutoff radius =",G_cutoff/au_length_aa," [1/A]"
    if(comm_is_root(nproc_id_global)) write(*,*) "# of cells used in R-cutoff=",ncell_max(:)
    if(comm_is_root(nproc_id_global)) write(*,*) "nc_max=",nc_max

    call init_dump_f
    call init_ewald_CRK

    do imol=1,nmol
       is = mol2species(imol)
       if(name_mol_s(is)=="WAT"   .or. &
          name_mol_s(is)=="WATER" .or. &
          name_mol_s(is)=="Water" .or. &
          name_mol_s(is)=="water" ) then
          is_wat = is
          exit
        endif
     enddo
     call get_Q0_K_wat( Q0ai,dQdr_wat,Kai,dKdr_wat,is_wat )
     Qai(:,:) = Q0ai(:,:)  !initial guess at the first time step

  end subroutine prepare_force_field_CRK


  subroutine cal_force_energy_CRK
    implicit none

    ! Reset forces and total energy
    Uene       = 0d0
    force(:,:) = 0d0

    call force_energy_intramolecular_CRK
    call force_energy_intermolecular_CRK
   !call force_energy_intermolecular_CRK_cluster

  end subroutine

!---------------------------------------------------------
  subroutine force_energy_intramolecular_CRK
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer :: ib_r, jb_r, ib_a, ib_rr, ib_ra, iatom1,iatom2,jatom1,jatom2, ipow, isign
    real(8) :: ene1, f_tmp(3), f_tmp_iatom1(3), f_tmp_iatom2(3), f_tmp_jatom1(3), f_tmp_jatom2(3)
    real(8) :: r12_vec(3,nb_r), r12_dist(nb_r), dr12_dist(nb_r), dr12_pow(1:6)
    real(8) :: r13_vec(3,nb_a), r13_dist(nb_a), dr13_dist(nb_a), dr12_c,dr12_c4,dr12_c5

    !calc dr for all 1-2 and 1-3 bonds
!$omp parallel do private(ib_r,iatom1,iatom2)
    do ib_r= 1,nb_r
       iatom1 = iatom1_b_r(ib_r)
       iatom2 = iatom2_b_r(ib_r)
       r12_vec(:,ib_r) = Rion(:,iatom1) - Rion(:,iatom2)
       r12_dist(ib_r)  = sqrt(sum(r12_vec(:,ib_r)*r12_vec(:,ib_r)))
       dr12_dist(ib_r) = r12_dist(ib_r) - r0_b_r(ib_r)
    enddo
!$omp end parallel do
!$omp parallel do private(ib_a,iatom1,iatom2)
    do ib_a= 1,nb_a
       iatom1 = iatom1_b_a(ib_a)
       iatom2 = iatom2_b_a(ib_a)
       r13_vec(:,ib_a) = Rion(:,iatom1) - Rion(:,iatom2)
       r13_dist(ib_a)  = sqrt(sum(r13_vec(:,ib_a)*r13_vec(:,ib_a)))
       dr13_dist(ib_a) = r13_dist(ib_a) - r0_b_a(ib_a)
    enddo
!$omp end parallel do

    ene1 = 0d0

    ! sum_n kn*(dr12)^n + kc*(dr12-rc)^5 term
!$omp parallel do private(ib_r,iatom1,iatom2,dr12_pow,ipow,isign,f_tmp) &
!$omp    private(dr12_c,dr12_c4,dr12_c5) &
!$omp    reduction(+:force,ene1)
    do ib_r= 1,nb_r
       iatom1 = iatom1_b_r(ib_r)
       iatom2 = iatom2_b_r(ib_r)
       dr12_pow(1)= dr12_dist(ib_r)

       do ipow=2,6
          dr12_pow(ipow)= dr12_pow(ipow-1) * dr12_pow(1)
       enddo

       if(dr12_pow(1).ge.0d0) then
          isign=1
       else
          isign=2
       endif

       do ipow=2,6
          f_tmp(:) = -ipow * dr12_pow(ipow-1) * r12_vec(:,ib_r)/r12_dist(ib_r)*kb_r(ib_r,ipow,isign)
          force(:,iatom1) = force(:,iatom1) + f_tmp(:)
          force(:,iatom2) = force(:,iatom2) - f_tmp(:)
          ene1 = ene1 + kb_r(ib_r,ipow,isign) * dr12_pow(ipow)
       enddo

       !this is just for preventing polarization catastrophy
       dr12_c = dr12_dist(ib_r) - r0_b_c(ib_r)
       if( dr12_c .gt. 0d0 ) then
          dr12_c4 = dr12_c**4d0
          dr12_c5 = dr12_c4*dr12_c
          f_tmp(:) = -5d0 * dr12_c4 * r12_vec(:,ib_r)/r12_dist(ib_r)*kb_c(ib_r)
          force(:,iatom1) = force(:,iatom1) + f_tmp(:)
          force(:,iatom2) = force(:,iatom2) - f_tmp(:)
          ene1 = ene1 + kb_c(ib_r) * dr12_c5
       endif
    enddo
!$omp end parallel do

    ! krr'*dr12*dr12' term
!$omp parallel do &
!$omp    private(ib_rr,ib_r,jb_r,iatom1,iatom2,jatom1,jatom2) &
!$omp    private(f_tmp_iatom1,f_tmp_iatom2,f_tmp_jatom1,f_tmp_jatom2) &
!$omp    reduction(+:force,ene1)
    do ib_rr= 1,nb_rr
       ib_r = ibond1_b_rr(ib_rr)
       jb_r = ibond2_b_rr(ib_rr)

       iatom1 = iatom1_b_r(ib_r)
       iatom2 = iatom2_b_r(ib_r)
       jatom1 = iatom1_b_r(jb_r)
       jatom2 = iatom2_b_r(jb_r)

       f_tmp_iatom1(:) = -r12_vec(:,ib_r)/r12_dist(ib_r) * dr12_dist(jb_r)*kb_rr(ib_rr)
       f_tmp_iatom2(:) = -f_tmp_iatom1(:)
       f_tmp_jatom1(:) = -r12_vec(:,jb_r)/r12_dist(jb_r) * dr12_dist(ib_r)*kb_rr(ib_rr)
       f_tmp_jatom2(:) = -f_tmp_jatom1(:)
       force(:,iatom1) = force(:,iatom1) + f_tmp_iatom1(:)
       force(:,iatom2) = force(:,iatom2) + f_tmp_iatom2(:)
       force(:,jatom1) = force(:,jatom1) + f_tmp_jatom1(:)
       force(:,jatom2) = force(:,jatom2) + f_tmp_jatom2(:)
       ene1 = ene1 + kb_rr(ib_rr) * dr12_dist(ib_r) * dr12_dist(jb_r)
    enddo
!$omp end parallel do

    ! (ktheta,r/2)*(dr13)^2 term
!$omp parallel do private(ib_a,iatom1,iatom2,f_tmp) reduction(+:force,ene1)
    do ib_a=1,nb_a
       iatom1 = iatom1_b_a(ib_a)
       iatom2 = iatom2_b_a(ib_a)
       f_tmp(:) = -kb_a(ib_a) * dr13_dist(ib_a) * r13_vec(:,ib_a)/r13_dist(ib_a)
       force(:,iatom1) = force(:,iatom1) + f_tmp(:)
       force(:,iatom2) = force(:,iatom2) - f_tmp(:)
       ene1 = ene1 + 0.5d0 * kb_a(ib_a) * dr13_dist(ib_a) * dr13_dist(ib_a)
    enddo
!$omp end parallel do

!$omp parallel do &
!$omp    private(ib_ra,ib_r,ib_a,iatom1,iatom2,jatom1,jatom2) &
!$omp    private(f_tmp_iatom1,f_tmp_iatom2,f_tmp_jatom1,f_tmp_jatom2) &
!$omp    reduction(+:force,ene1)
    do ib_ra=1,nb_ra
       ib_r = ibond_b_ra(ib_ra)
       ib_a = iangle_b_ra(ib_ra)

       iatom1 = iatom1_b_r(ib_r)
       iatom2 = iatom2_b_r(ib_r)
       jatom1 = iatom1_b_a(ib_a)
       jatom2 = iatom2_b_a(ib_a)

       f_tmp_iatom1(:) = -r12_vec(:,ib_r)/r12_dist(ib_r) * dr13_dist(ib_a)*kb_ra(ib_ra)
       f_tmp_iatom2(:) = -f_tmp_iatom1(:)
       f_tmp_jatom1(:) = -r13_vec(:,ib_a)/r13_dist(ib_a) * dr12_dist(ib_r)*kb_ra(ib_ra)
       f_tmp_jatom2(:) = -f_tmp_jatom1(:)
       force(:,iatom1) = force(:,iatom1) + f_tmp_iatom1(:)
       force(:,iatom2) = force(:,iatom2) + f_tmp_iatom2(:)
       force(:,jatom1) = force(:,jatom1) + f_tmp_jatom1(:)
       force(:,jatom2) = force(:,jatom2) + f_tmp_jatom2(:)
       ene1 = ene1 + kb_ra(ib_ra) * dr12_dist(ib_r) * dr13_dist(ib_a)
    enddo
!$omp end parallel do

    Uene = Uene + ene1

  end subroutine

!---------------------------------------------------------
  subroutine force_energy_intermolecular_CRK
    use salmon_math
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer :: ia,ia_wX, ja_wX, imol,jmol, is,is_wat,index_image(3,nc_max),nc,ic
    integer :: iatom_c, jatom_c, iatom, iter, iX, ixyz, iG
    real(8) :: drij(3),drij0(3),rij,rij_inv,rij2,rij2_inv,rij6_inv,rij12_inv,r2_cutoff
    real(8) :: ene_vdw, ene_vdw6, ene_vdw12
    real(8) :: eps_sgm6_vdw_wat, eps_sgm12_vdw_wat
    real(8) :: ene_cln, ene_cln_real, ene_cln_reci, ene_cln_self, ene_cln_mask
    real(8) :: e_tmp, f_tmp, f_tmp_iatom(3), dump_f, dump_dfdr(3)
    real(8) :: force_wX(3,4,nmol), dRXdRatom(3,2,3,3,nmol)
    real(8) :: Qai_old(4,nmol),Q0ai(4,nmol), dQave
    real(8) :: QiQj, dQdr_wat(3,4,3,nmol)
    real(8) :: Kai(4,4,nmol), dKdr_wat(3,4,4,3,nmol) 
    real(8) :: Vai(4,nmol), tmp_sumKV(3,4)
    real(8) :: erf_rij,exp_rij, a_ewald_rij
    real(8) :: tmp_s,tmp_c,tmp_coef,tmp1,tmp2,tmp3

    if(nmol_s.ge.2) then
       call Err_finalize("Now combination rule of vdw is not supported except for water")
       stop
    endif

    !for water (oxygen)
    do imol=1,nmol
       is = mol2species(imol)
       if(name_mol_s(is)=="WAT"   .or. &
          name_mol_s(is)=="WATER" .or. &
          name_mol_s(is)=="Water" .or. &
          name_mol_s(is)=="water" ) then
          is_wat = is
          ia = mol2atom_top(imol)
          eps_sgm6_vdw_wat = eps_vdw(ia) * sgm_vdw(ia)**6d0
          eps_sgm12_vdw_wat= eps_vdw(ia) * sgm_vdw(ia)**12d0
          exit
        endif
    enddo

    r2_cutoff = r_cutoff * r_cutoff

    !generate X sites for each water molecules
    call get_Xatom_wat(Rxyz_wX,dRXdRatom)
    
    !get Q0 and Kai
    call get_Q0_K_wat( Q0ai,dQdr_wat,Kai,dKdr_wat,is_wat )
    
    !calculate and save some terms for Ewald
    ene_vdw6  = 0d0
    ene_vdw12 = 0d0
    ene_cln   = 0d0

    dump_f       = 1d0
    dump_dfdr(:) = 0d0

    !Qai(:,:) = Q0ai(:,:)
    Vai(:,:) = 0d0

!$omp parallel do collapse(2) &
!$omp    private(imol,jmol,iatom_c,jatom_c,ic,nc) & 
!$omp    private(drij,drij0,index_image,rij2,rij2_inv,rij6_inv,rij12_inv) &
!$omp    private(f_tmp,f_tmp_iatom,ia_wX,ja_wX,rij,rij_inv) &
!$omp    private(a_ewald_rij,erf_rij,dump_f,dump_dfdr,tmp1) &
!$omp    reduction(+:Vai,ene_vdw6,ene_vdw12,force,force_wX)
    do imol=1,nmol
    do jmol=1,nmol
       iatom_c = mol2atom_cnt(imol)
       jatom_c = mol2atom_cnt(jmol)

       drij0(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
       call image_periodic(drij0,index_image,nc)
       do ic=1,nc
          drij(:) = drij0(:) + index_image(:,ic)*aL(:)
          rij2 = sum(drij(:)*drij(:))
          if(rij2 .gt. r2_cutoff) cycle   !cut-off

          if(imol==jmol) then
             if(rij2.le.1d-5) cycle  !exclude same molecule in the same cell
          endif

          !! VDW interaction (water, only between oxygen atoms)
          rij2_inv = 1d0/rij2
          rij6_inv = rij2_inv*rij2_inv*rij2_inv
          rij12_inv= rij6_inv*rij6_inv
          ene_vdw6  = ene_vdw6  + rij6_inv   *0.5d0
          ene_vdw12 = ene_vdw12 + rij12_inv  *0.5d0
          f_tmp = 2d0*eps_sgm12_vdw_wat*rij12_inv - eps_sgm6_vdw_wat*rij6_inv
          f_tmp_iatom(:) = -f_tmp * 24d0 * rij2_inv * drij(:) *0.5d0
          force(:,iatom_c) = force(:,iatom_c) + f_tmp_iatom(:)
          force(:,jatom_c) = force(:,jatom_c) - f_tmp_iatom(:)

          !! Coulomb interaction  Here, just get Vabi before iteration
          !!    Real Space Term of Ewald (the other terms are after the loop)
          do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          do ja_wX= 1,natom_wX_mol(jmol)
             drij(:) = (Rxyz_wX(:,ja_wX,jmol)-Rxyz_wX(:,ia_wX,imol)) + index_image(:,ic)*aL(:)
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij

             a_ewald_rij = a_ewald * rij
             erf_rij     = erf_salmon(a_ewald_rij)

             call get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,0)
             tmp1 = ( dump_f - erf_rij )*rij_inv  *0.5d0
             Vai(ia_wX,imol) = Vai(ia_wX,imol) + Qai(ja_wX,jmol) * tmp1
             Vai(ja_wX,jmol) = Vai(ja_wX,jmol) + Qai(ia_wX,imol) * tmp1

          enddo
          enddo
       enddo !ic

    enddo   !jmol
    enddo   !imol
!$omp end parallel do


    !! Reciprocal Space Term of Ewald
    call get_sumQsinGr_sumQcosGr_ewald(Rxyz_wX,Qai)
    tmp1 = 1d0/(4d0*a_ewald*a_ewald)
    tmp2 = 2d0*pi/Vuc * 2d0
!$omp parallel do &
!$omp private(iG,tmp_coef,imol,ia_wX,tmp3,tmp_c,tmp_s) reduction(+:Vai)
    do iG=1,nG_cmd
       tmp_coef = tmp2 * exp(-G2(iG)*tmp1)/G2(iG)
       do imol=1,nmol
       do ia_wX=1,natom_wX_mol(imol)
          tmp3  = sum(Gvec(:,iG)*Rxyz_wX(:,ia_wX,imol))
          tmp_c = cos(tmp3) * sumQicosGri(iG)
          tmp_s = sin(tmp3) * sumQisinGri(iG)
          Vai(ia_wX,imol) = Vai(ia_wX,imol) + tmp_coef *(tmp_c + tmp_s)
       end do
       end do
    end do
!$omp end parallel do

    !! Self Term of Ewald
    !(for the same atom)
    tmp1 = a_ewald/sqrt(pi) * 2d0
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:Vai)
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       Vai(ia_wX,imol) = Vai(ia_wX,imol) - tmp1 * Qai(ia_wX,imol)
    enddo
    enddo
!$omp end parallel do

    !(for the different atoms in the same molecule)
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX,ja_wX,drij,rij2,rij,rij_inv,erf_rij) reduction(+:Vai)
    do imol=1,nmol
       do ia_wX= 1,natom_wX_mol(imol)
       do ja_wX= ia_wX+1,natom_wX_mol(imol)

          drij(:) = Rxyz_wX(:,ja_wX,imol) - Rxyz_wX(:,ia_wX,imol)
          rij2    = sum(drij(:)*drij(:))
          rij     = sqrt(rij2)
          rij_inv = 1d0/rij
          erf_rij = erf_salmon(a_ewald*rij)

          Vai(ia_wX,imol) = Vai(ia_wX,imol) - Qai(ja_wX,imol) * erf_rij * rij_inv
          Vai(ja_wX,imol) = Vai(ja_wX,imol) - Qai(ia_wX,imol) * erf_rij * rij_inv

       enddo
       enddo
    enddo
!$omp end parallel do

    !Add external electric field
    if(flag_use_Et)then
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:Vai)
    do imol=1,nmol
    do ia_wX= 1,natom_wX_mol(imol)
       Vai(ia_wX,imol) = Vai(ia_wX,imol) - sum(Rxyz_wX(:,ia_wX,imol)*Et(:))
    enddo
    enddo
!$omp end parallel do
    endif

    !! Iteration for Coulomb
    do iter=1,iter_max

       !get charge Qai using Vai
       Qai_old(:,:) = Qai(:,:)
       dQave = 0d0
       do imol=1,nmol
       do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          Qai(ia_wX,imol) = Q0ai(ia_wX,imol) + sum( Kai(ia_wX,:,imol)*Vai(:,imol) )
          dQave = dQave + ( Qai(ia_wX,imol) - Qai_old(ia_wX,imol) )**2
       enddo
       enddo
       dQave = sqrt( dQave/NI_wX )

       !get Vai using Qai
       Vai(:,:) = 0d0

!$omp parallel do collapse(2) &
!$omp    private(imol,jmol,iatom_c,jatom_c,ic,nc) & 
!$omp    private(drij,drij0,index_image,rij2,rij2_inv,rij,rij_inv) &
!$omp    private(ia_wX,ja_wX,a_ewald_rij,erf_rij,dump_f,dump_dfdr,tmp1) &
!$omp    reduction(+:Vai)
       do imol=1,nmol
       do jmol=1,nmol
          iatom_c = mol2atom_cnt(imol)
          jatom_c = mol2atom_cnt(jmol)

          drij0(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
          call image_periodic(drij0,index_image,nc)
          do ic=1,nc
          drij(:) = drij0(:) + index_image(:,ic)*aL(:)
          rij2 = sum(drij(:)*drij(:))
          if(rij2 .gt. r2_cutoff) cycle   !cut-off

          if(imol==jmol) then
             if(rij2.le.1d-5) cycle  !exclude same molecule in the same cell
          endif

          !! Coulomb interaction  Here, just get Vabi before iteration
          !!    Real Space Term of Ewald (the other terms are after the loop)
          do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          do ja_wX= 1,natom_wX_mol(jmol)
             drij(:) = (Rxyz_wX(:,ja_wX,jmol)-Rxyz_wX(:,ia_wX,imol)) + index_image(:,ic)*aL(:)
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij

             a_ewald_rij = a_ewald * rij
             erf_rij     = erf_salmon(a_ewald_rij)

             call get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,0)
             tmp1 = ( dump_f - erf_rij )*rij_inv *0.5d0
             Vai(ia_wX,imol) = Vai(ia_wX,imol) + Qai(ja_wX,jmol) * tmp1
             Vai(ja_wX,jmol) = Vai(ja_wX,jmol) + Qai(ia_wX,imol) * tmp1
          enddo
          enddo
          enddo  !ic
       enddo  !jmol
       enddo  !imol
!$omp end parallel do

       !! Reciprocal Space Term of Ewald
       call get_sumQsinGr_sumQcosGr_ewald(Rxyz_wX,Qai)
       tmp1 = 1d0/(4d0*a_ewald*a_ewald)
       tmp2 = 2d0*pi/Vuc * 2d0
!$omp parallel do &
!$omp private(iG,tmp_coef,imol,ia_wX,tmp3,tmp_c,tmp_s) reduction(+:Vai)
       do iG=1,nG_cmd
          tmp_coef = tmp2 * exp(-G2(iG)*tmp1)/G2(iG)
          do imol=1,nmol
          do ia_wX=1,natom_wX_mol(imol)
             tmp3  = sum(Gvec(:,iG)*Rxyz_wX(:,ia_wX,imol))
             tmp_c = cos(tmp3) * sumQicosGri(iG)
             tmp_s = sin(tmp3) * sumQisinGri(iG)
             Vai(ia_wX,imol) = Vai(ia_wX,imol) + tmp_coef *(tmp_c + tmp_s)
          end do
          end do
       end do
!$omp end parallel do

       !! Self Term of Ewald
       !(for the same atom)
       tmp1 = -a_ewald/sqrt(pi) * 2d0
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:Vai)
       do imol=1,nmol
       do ia_wX=1,natom_wX_mol(imol)
          Vai(ia_wX,imol) = Vai(ia_wX,imol) + tmp1 * Qai(ia_wX,imol)
       enddo
       enddo
!$omp end parallel do

       !(for the different atoms in the same molecule)
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX,ja_wX,drij,rij2,rij,rij_inv,erf_rij) reduction(+:Vai)
       do imol=1,nmol
          do ia_wX= 1,natom_wX_mol(imol)
          do ja_wX= ia_wX+1,natom_wX_mol(imol)

             drij(:) = Rxyz_wX(:,ja_wX,imol) - Rxyz_wX(:,ia_wX,imol)
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij
             erf_rij = erf_salmon(a_ewald*rij)

             Vai(ia_wX,imol) = Vai(ia_wX,imol) - Qai(ja_wX,imol) * erf_rij * rij_inv
             Vai(ja_wX,imol) = Vai(ja_wX,imol) - Qai(ia_wX,imol) * erf_rij * rij_inv

          enddo
          enddo
       enddo
!$omp end parallel do

       !Add external electric field
       if(flag_use_Et)then
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:Vai)
       do imol=1,nmol
       do ia_wX= 1,natom_wX_mol(imol)
          Vai(ia_wX,imol) = Vai(ia_wX,imol) - sum(Rxyz_wX(:,ia_wX,imol)*Et(:))
       enddo
       enddo
!$omp end parallel do
       endif

       !check convergence
       if(dQave .le. dQave_thresh) then
         !write(*,*) "#Qai is converged: dQave=", dQave, iter
          exit
       else if(iter == iter_max) then
          if (comm_is_root(nproc_id_global)) &
          write(*,*) "#iter=iter_max: not converged", dQave, iter
       endif

    enddo  !iter

    !log if necessary: test hoge
    !write(*,*) "#abs(Qmax)=", maxval(abs(Qai(:,:)))
    !write(*,*) "#ave of Qai@Xsite=",(sum(Qai(1,:))+sum(Qai(2,:)))/(2*nmol)
    !write(*,*) "#ave of Qai@Hatom=",(sum(Qai(3,:))+sum(Qai(4,:)))/(2*nmol)


    ! Calculate Energy and Focer of Coulomb interaction after convergence
    force_wX(:,:,:) = 0d0
    ene_cln_real = 0d0
    ene_cln_reci = 0d0
    ene_cln_self = 0d0
    ene_cln_mask = 0d0
!$omp parallel do collapse(2) &
!$omp    private(imol,jmol,iatom_c,jatom_c,ia_wX,ja_wX,ic,nc) & 
!$omp    private(drij,drij0,index_image,rij2,rij2_inv,rij,rij_inv) &
!$omp    private(QiQj,a_ewald_rij,erf_rij,exp_rij,dump_f,dump_dfdr,tmp1) &
!$omp    private(e_tmp,f_tmp_iatom,ia,ixyz,tmp_sumKV,iatom) &
!$omp    reduction(+:force_wX,force,ene_cln_real)
    do imol=1,nmol
    do jmol=1,nmol
       iatom_c = mol2atom_cnt(imol)
       jatom_c = mol2atom_cnt(jmol)

       drij0(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
       call image_periodic(drij0,index_image,nc)
       do ic=1,nc
          drij(:) = drij0(:) + index_image(:,ic)*aL(:)
          rij2 = sum(drij(:)*drij(:))
          if(rij2 .gt. r2_cutoff) cycle   !cut-off

          if(imol==jmol) then
             if(rij2.le.1d-5) cycle  !exclude same molecule in the same cell
          endif

          !! Coulomb interaction 
          !!    Real Space Term of Ewald (the other terms are after the loop)
          do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          do ja_wX= 1,natom_wX_mol(jmol)
             drij(:) = (Rxyz_wX(:,ja_wX,jmol)-Rxyz_wX(:,ia_wX,imol)) + index_image(:,ic)*aL(:)
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij
             rij2_inv= 1d0/rij2

             QiQj        = Qai(ia_wX,imol) * Qai(ja_wX,jmol)
             a_ewald_rij = a_ewald * rij
             erf_rij    = erf_salmon(a_ewald_rij)
             exp_rij     = exp(-a_ewald_rij * a_ewald_rij)

             call get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,1)
             e_tmp          = QiQj * ( dump_f -  erf_rij )*rij_inv
             f_tmp_iatom(:) = QiQj *( -rij_inv*( drij(:) * rij2_inv * dump_f + dump_dfdr(:) ) &
                                      +drij(:)*rij2_inv * (erf_rij * rij_inv - coef_derfc*exp_rij) )
             e_tmp          = e_tmp * 0.5d0
             f_tmp_iatom(:) = f_tmp_iatom(:) * 0.5d0

             force_wX(:,ia_wX,imol) = force_wX(:,ia_wX,imol) + f_tmp_iatom(:)
             force_wX(:,ja_wX,jmol) = force_wX(:,ja_wX,jmol) - f_tmp_iatom(:)
             ene_cln_real = ene_cln_real + e_tmp
          enddo
          enddo
       enddo    !ic
    enddo    ! jmol
    enddo    ! imol
!$omp end parallel do

!$omp parallel do &
!$omp    private(imol,ia,ia_wX,ixyz,e_tmp,f_tmp_iatom,tmp_sumKV,iatom) &
!$omp    reduction(+:force_wX,force,ene_cln)
    do imol=1,nmol
       !! Add force coming from conformation-dependent partial charge
       do ia= 1,natom_mol(imol)   !O,H1,H2
          do ixyz= 1,3
             do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
                tmp_sumKV(ixyz,ia_wX) = sum( dKdr_wat(ixyz,ia_wX,:,ia,imol) * Vai(:,imol) )
             enddo
             f_tmp_iatom(ixyz) = -0.5d0* sum(tmp_sumKV(ixyz,:)*Vai(:,imol)) - sum(dQdr_wat(ixyz,:,ia,imol)*Vai(:,imol))
          enddo
          iatom = mol2atom_top(imol) + ia-1
          force(:,iatom) = force(:,iatom) + f_tmp_iatom(:)
       enddo

       !! Add energy coming from sum(KabVaiVbi) term
       e_tmp = 0d0
       do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          e_tmp = e_tmp + sum(Kai(ia_wX,:,imol)*Vai(:,imol)) * Vai(ia_wX,imol)
       enddo
       ene_cln = ene_cln - 0.5d0*e_tmp
    enddo      ! imol
!$omp end parallel do


    !! Reciprocal Space Term of Ewald
    !call get_sumQsinGr_sumQcosGr_ewald(Rxyz_wX,Qai)
    tmp1 = 1d0/(4d0*a_ewald*a_ewald)
    tmp2 = 2d0*pi/Vuc
!$omp parallel do &
!$omp   private(iG,tmp_coef,e_tmp,imol,ia_wX,tmp3,tmp_c,tmp_s,f_tmp_iatom) &
!$omp   reduction(+:force_wX,ene_cln_reci)
    do iG=1,nG_cmd
       tmp_coef = tmp2 * exp(-G2(iG)*tmp1)/G2(iG)
       e_tmp = tmp_coef * (sumQicosGri(iG)**2 + sumQisinGri(iG)**2)
       ene_cln_reci = ene_cln_reci + e_tmp
       do imol=1,nmol
       do ia_wX=1,natom_wX_mol(imol)
          tmp3  = sum(Gvec(:,iG)*Rxyz_wX(:,ia_wX,imol))
          tmp_c = Qai(ia_wX,imol) * cos(tmp3)
          tmp_s = Qai(ia_wX,imol) * sin(tmp3)
          f_tmp_iatom(:) = 2d0 * tmp_coef * Gvec(:,iG) &
                  * ( tmp_s * sumQicosGri(iG) - tmp_c * sumQisinGri(iG) )
          force_wX(:,ia_wX,imol) = force_wX(:,ia_wX,imol) + f_tmp_iatom(:)
       end do
       end do
    end do
!$omp end parallel do


    !! Self Term of Ewald
    !(for the same atom)
    ene_cln_self = 0d0
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:ene_cln_self)
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       ene_cln_self = ene_cln_self + Qai(ia_wX,imol) * Qai(ia_wX,imol)
    enddo
    enddo
!$omp end parallel do
    ene_cln_self = -ene_cln_self * a_ewald/sqrt(pi)


    !(for the different atoms in the same molecule)
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX,ja_wX,drij,rij2,rij2_inv,rij,rij_inv) &
!$omp private(QiQj,erf_rij,exp_rij,e_tmp,f_tmp_iatom) &
!$omp reduction(+:force_wX,ene_cln_mask)
    do imol=1,nmol
       do ia_wX= 1,natom_wX_mol(imol)
       do ja_wX= ia_wX+1,natom_wX_mol(imol)

          drij(:) = Rxyz_wX(:,ja_wX,imol) - Rxyz_wX(:,ia_wX,imol)
          rij2    = sum(drij(:)*drij(:))
          rij     = sqrt(rij2)
          rij2_inv= 1d0/rij2
          rij_inv = 1d0/rij

          QiQj    = Qai(ia_wX,imol) * Qai(ja_wX,imol)
          erf_rij = erf_salmon(a_ewald*rij)
          exp_rij = exp(-a_ewald*a_ewald*rij2)

          e_tmp = -QiQj * erf_rij * rij_inv
          f_tmp_iatom(:) = drij(:)*rij2_inv * QiQj *(erf_rij * rij_inv - coef_derfc*exp_rij) 

          force_wX(:,ia_wX,imol) = force_wX(:,ia_wX,imol) + f_tmp_iatom(:)
          force_wX(:,ja_wX,imol) = force_wX(:,ja_wX,imol) - f_tmp_iatom(:)
          ene_cln_mask = ene_cln_mask + e_tmp
       enddo
       enddo
    enddo
!$omp end parallel do

    ! Add force by electric field
    if(flag_use_Et)then
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       force_wX(:,ia_wX,imol) = force_wX(:,ia_wX,imol) + Qai(ia_wX,imol) * Et(:)
    enddo
    enddo
    endif

    ! convert from system with X atoms into ordinary real atom system
    do imol=1,nmol
       !force on H1,H2 (ordinary Coulomb interaction)
       do ia= 2,3  !H1,H2
          iatom = mol2atom_top(imol) + ia-1
          ia_wX = ia+1
          force(:,iatom) = force(:,iatom) + force_wX(:,ia_wX,imol)
       enddo

       ! force through X1,X2 sites in imol & jmol
       do ia= 1,3    ! O,H1,H2 atoms in imol
          iatom = mol2atom_top(imol) + ia-1
          do iX= 1,2       ! X1,X2 sites in imol/jmol
          do ixyz=1,3
             force(ixyz,iatom) = force(ixyz,iatom) + sum(dRXdRatom(:,iX,ixyz,ia,imol)*force_wX(:,iX,imol))
          enddo
          enddo
       enddo
    enddo     ! imol

    ene_vdw = 4d0*(eps_sgm12_vdw_wat * ene_vdw12 - eps_sgm6_vdw_wat * ene_vdw6)
    ene_cln = ene_cln + ene_cln_real + ene_cln_reci + ene_cln_self + ene_cln_mask
    Uene = Uene + ene_cln + ene_vdw

    !!check of coulomb energy from different calculation
    !e_tmp=0d0
    !do imol=1,nmol
    !do ia_wX=1,natom_wX_mol(imol)
    !   e_tmp = e_tmp + 0.5d0*Qai(ia_wX,imol)*Vai(ia_wX,imol)
    !enddo
    !enddo
    !do imol=1,nmol
    !do ia_wX=1,natom_wX_mol(imol)
    !do ja_wX=1,natom_wX_mol(imol)
    !   e_tmp = e_tmp - 0.5d0*Kai(ia_wX,ja_wX,imol)*Vai(ia_wX,imol)*Vai(ja_wX,imol)
    !enddo
    !enddo
    !enddo
    !write(*,'(a,3e22.10)') "#check-coulomb-energy:", e_tmp, ene_cln, e_tmp-ene_cln

!hoge
!    !write for check 
!    write(*,*) "------ output for check -------"
!    do imol=1,nmol
!    do ia=1,3
!       iatom = mol2atom_top(imol) + ia-1
!       write(*,5) imol,iatom,(Rion(ixyz,iatom),ixyz=1,3)
!    enddo
!    enddo
!  5 format("Rion=",2i4,3e20.10)
!    write(*,10) Uene, (Uene-ene_cln-ene_vdw), ene_vdw, ene_cln
! 10 format("Energy=", 4e20.10)
!    do imol=1,nmol
!       write(*,20) imol,(Q0ai(ia_wX,imol),ia_wX=1,4)
!    enddo
! 20 format("Charge(Q0)=",i4,4e20.10)
!    do imol=1,nmol
!       write(*,25) imol,(Qai(ia_wX,imol),ia_wX=1,4)
!    enddo
! 25 format("Charge(Q)=",i4,4e20.10)
!    do imol=1,nmol
!    do ia=1,3
!       iatom = mol2atom_top(imol) + ia-1
!       write(*,30) imol,iatom,(force(ixyz,iatom),ixyz=1,3)
!    enddo
!    enddo
! 30 format("Force=",2i4,3e20.10)
!    write(*,*) "-------------------------------"
!    stop

  end subroutine force_energy_intermolecular_CRK

!---------------------------------------------------------
  subroutine force_energy_intermolecular_CRK_cluster
    use salmon_math
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer :: ia,ia_wX, ja_wX, imol,jmol, is,is_wat
    integer :: iatom_c, jatom_c, iatom, iter, iX, ixyz
    real(8) :: drij(3),rij,rij_inv,rij2,rij2_inv,rij6_inv,rij12_inv
    real(8) :: ene_vdw, ene_vdw6, ene_vdw12
    real(8) :: ene_cln, eps_sgm6_vdw_wat, eps_sgm12_vdw_wat
    real(8) :: e_tmp, f_tmp, f_tmp_iatom(3), dump_f, dump_dfdr(3)
    real(8) :: force_wX(3,4,nmol), dRXdRatom(3,2,3,3,nmol)
    real(8) :: Qai_old(4,nmol),Q0ai(4,nmol), dQave
    real(8) :: QiQj, dQdr_wat(3,4,3,nmol)
    real(8) :: Kai(4,4,nmol), dKdr_wat(3,4,4,3,nmol) 
    real(8) :: Vai(4,nmol), tmp_sumKV(3,4), tmp1


    if(nmol_s.ge.2) then
       call Err_finalize("Now combination rule of vdw is not supported except for water")
       stop
    endif

    !for water (oxygen)
    do imol=1,nmol
       is = mol2species(imol)
       if(name_mol_s(is)=="WAT"   .or. &
          name_mol_s(is)=="WATER" .or. &
          name_mol_s(is)=="Water" .or. &
          name_mol_s(is)=="water" ) then
          is_wat = is
          ia = mol2atom_top(imol)
          eps_sgm6_vdw_wat = eps_vdw(ia) * sgm_vdw(ia)**6d0
          eps_sgm12_vdw_wat= eps_vdw(ia) * sgm_vdw(ia)**12d0
          exit
        endif
    enddo

    !generate X sites for each water molecules
    call get_Xatom_wat(Rxyz_wX,dRXdRatom)
    
    !get Q0 and Kai
    call get_Q0_K_wat( Q0ai,dQdr_wat,Kai,dKdr_wat,is_wat )
    
    !calculate and save some terms for Ewald
    ene_vdw6  = 0d0
    ene_vdw12 = 0d0
    ene_cln   = 0d0

    dump_f       = 1d0
    dump_dfdr(:) = 0d0

    !Qai(:,:) = Q0ai(:,:)
    Vai(:,:) = 0d0

!!$omp parallel do collapse(2) &
!!$omp    private(imol,jmol,iatom_c,jatom_c,ic,nc) & 
!!$omp    private(drij,drij0,index_image,rij2,rij2_inv,rij6_inv,rij12_inv) &
!!$omp    private(f_tmp,f_tmp_iatom,ia_wX,ja_wX,rij,rij_inv) &
!!$omp    private(a_ewald_rij,erf_rij,dump_f,dump_dfdr,tmp1) &
!!$omp    reduction(+:Vai,ene_vdw6,ene_vdw12,force,force_wX)
    do imol=1,nmol
    do jmol=1,nmol

       if(imol==jmol) cycle
       iatom_c = mol2atom_cnt(imol)
       jatom_c = mol2atom_cnt(jmol)

       drij(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
       rij2 = sum(drij(:)*drij(:))

          !! VDW interaction (water, only between oxygen atoms)
          rij2_inv = 1d0/rij2
          rij6_inv = rij2_inv*rij2_inv*rij2_inv
          rij12_inv= rij6_inv*rij6_inv
          ene_vdw6  = ene_vdw6  + rij6_inv   *0.5d0
          ene_vdw12 = ene_vdw12 + rij12_inv  *0.5d0
          f_tmp = 2d0*eps_sgm12_vdw_wat*rij12_inv - eps_sgm6_vdw_wat*rij6_inv
          f_tmp_iatom(:) = -f_tmp * 24d0 * rij2_inv * drij(:) *0.5d0
          force(:,iatom_c) = force(:,iatom_c) + f_tmp_iatom(:)
          force(:,jatom_c) = force(:,jatom_c) - f_tmp_iatom(:)

          !! Coulomb interaction  Here, just get Vabi before iteration
          !!    Real Space Term of Ewald (the other terms are after the loop)
          do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          do ja_wX= 1,natom_wX_mol(jmol)
             drij(:) = (Rxyz_wX(:,ja_wX,jmol)-Rxyz_wX(:,ia_wX,imol))
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij

             call get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,0)
             tmp1 = dump_f * rij_inv  * 0.5d0
             Vai(ia_wX,imol) = Vai(ia_wX,imol) + Qai(ja_wX,jmol) * tmp1
             Vai(ja_wX,jmol) = Vai(ja_wX,jmol) + Qai(ia_wX,imol) * tmp1

          enddo
          enddo

    enddo   !jmol
    enddo   !imol
!!$omp end parallel do


    !Add external electric field
    if(flag_use_Et)then
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:Vai)
    do imol=1,nmol
    do ia_wX= 1,natom_wX_mol(imol)
       Vai(ia_wX,imol) = Vai(ia_wX,imol) - sum(Rxyz_wX(:,ia_wX,imol)*Et(:))
    enddo
    enddo
!$omp end parallel do
    endif

    !! Iteration for Coulomb
    do iter=1,iter_max

       !get charge Qai using Vai
       Qai_old(:,:) = Qai(:,:)
       dQave = 0d0
       do imol=1,nmol
       do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          Qai(ia_wX,imol) = Q0ai(ia_wX,imol) + sum( Kai(ia_wX,:,imol)*Vai(:,imol) )
          dQave = dQave + ( Qai(ia_wX,imol) - Qai_old(ia_wX,imol) )**2
       enddo
       enddo
       dQave = sqrt( dQave/NI_wX )

       !get Vai using Qai
       Vai(:,:) = 0d0

!!$omp parallel do collapse(2) &
!!$omp    private(imol,jmol,iatom_c,jatom_c,ic,nc) & 
!!$omp    private(drij,drij0,index_image,rij2,rij2_inv,rij,rij_inv) &
!!$omp    private(ia_wX,ja_wX,a_ewald_rij,erf_rij,dump_f,dump_dfdr,tmp1) &
!!$omp    reduction(+:Vai)
       do imol=1,nmol
       do jmol=1,nmol

          if(imol==jmol) cycle
         !iatom_c = mol2atom_cnt(imol)
         !jatom_c = mol2atom_cnt(jmol)

         !drij(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
         !rij2 = sum(drij(:)*drij(:))

          !! Coulomb interaction  Here, just get Vabi before iteration
          !!    Real Space Term of Ewald (the other terms are after the loop)
          do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          do ja_wX= 1,natom_wX_mol(jmol)
             drij(:) = (Rxyz_wX(:,ja_wX,jmol)-Rxyz_wX(:,ia_wX,imol))
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij

             call get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,0)
             tmp1 = dump_f *rij_inv *0.5d0
             Vai(ia_wX,imol) = Vai(ia_wX,imol) + Qai(ja_wX,jmol) * tmp1
             Vai(ja_wX,jmol) = Vai(ja_wX,jmol) + Qai(ia_wX,imol) * tmp1
          enddo
          enddo
       enddo  !jmol
       enddo  !imol
!!$omp end parallel do

       !Add external electric field
       if(flag_use_Et)then
!$omp parallel do collapse(2) &
!$omp private(imol,ia_wX) reduction(+:Vai)
       do imol=1,nmol
       do ia_wX= 1,natom_wX_mol(imol)
          Vai(ia_wX,imol) = Vai(ia_wX,imol) - sum(Rxyz_wX(:,ia_wX,imol)*Et(:))
       enddo
       enddo
!$omp end parallel do
       endif

       !check convergence
       if(dQave .le. dQave_thresh) then
         !write(*,*) "#Qai is converged: dQave=", dQave, iter
          exit
       else if(iter == iter_max) then
          if (comm_is_root(nproc_id_global)) &
          write(*,*) "#iter=iter_max: not converged", dQave, iter
       endif

    enddo  !iter

    !log if necessary: test hoge
    !write(*,*) "#abs(Qmax)=", maxval(abs(Qai(:,:)))
    !write(*,*) "#ave of Qai@Xsite=",(sum(Qai(1,:))+sum(Qai(2,:)))/(2*nmol)
    !write(*,*) "#ave of Qai@Hatom=",(sum(Qai(3,:))+sum(Qai(4,:)))/(2*nmol)


    ! Calculate Energy and Focer of Coulomb interaction after convergence
    force_wX(:,:,:) = 0d0
!!$omp parallel do collapse(2) &
!!$omp    private(imol,jmol,iatom_c,jatom_c,ia_wX,ja_wX,ic,nc) & 
!!$omp    private(drij,drij0,index_image,rij2,rij2_inv,rij,rij_inv) &
!!$omp    private(QiQj,a_ewald_rij,erf_rij,exp_rij,dump_f,dump_dfdr,tmp1) &
!!$omp    private(e_tmp,f_tmp_iatom,ia,ixyz,tmp_sumKV,iatom) &
!!$omp    reduction(+:force_wX,force,ene_cln_real)
    do imol=1,nmol
    do jmol=1,nmol

       if(imol==jmol) cycle
       !iatom_c = mol2atom_cnt(imol)
       !jatom_c = mol2atom_cnt(jmol)

       !drij(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
       !rij2 = sum(drij(:)*drij(:))

          !! Coulomb interaction 
          !!    Real Space Term of Ewald (the other terms are after the loop)
          do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          do ja_wX= 1,natom_wX_mol(jmol)
             drij(:) = (Rxyz_wX(:,ja_wX,jmol)-Rxyz_wX(:,ia_wX,imol))
             rij2    = sum(drij(:)*drij(:))
             rij     = sqrt(rij2)
             rij_inv = 1d0/rij
             rij2_inv= 1d0/rij2

             QiQj    = Qai(ia_wX,imol) * Qai(ja_wX,jmol)

             call get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,1)
             e_tmp          =  QiQj * rij_inv * dump_f
             f_tmp_iatom(:) = -QiQj * rij_inv*( drij(:) * rij2_inv * dump_f + dump_dfdr(:) )
             e_tmp          = e_tmp * 0.5d0
             f_tmp_iatom(:) = f_tmp_iatom(:) * 0.5d0

             force_wX(:,ia_wX,imol) = force_wX(:,ia_wX,imol) + f_tmp_iatom(:)
             force_wX(:,ja_wX,jmol) = force_wX(:,ja_wX,jmol) - f_tmp_iatom(:)
             ene_cln = ene_cln + e_tmp
          enddo
          enddo
    enddo    ! jmol
    enddo    ! imol
!!$omp end parallel do


!$omp parallel do &
!$omp    private(imol,ia,ia_wX,ixyz,e_tmp,f_tmp_iatom,tmp_sumKV,iatom) &
!$omp    reduction(+:force_wX,force,ene_cln)
    do imol=1,nmol
       !! Add force coming from conformation-dependent partial charge
       do ia= 1,natom_mol(imol)   !O,H1,H2
          do ixyz= 1,3
             do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
                tmp_sumKV(ixyz,ia_wX) = sum( dKdr_wat(ixyz,ia_wX,:,ia,imol) * Vai(:,imol) )
             enddo
             f_tmp_iatom(ixyz) = -0.5d0* sum(tmp_sumKV(ixyz,:)*Vai(:,imol)) - sum(dQdr_wat(ixyz,:,ia,imol)*Vai(:,imol))
          enddo
          iatom = mol2atom_top(imol) + ia-1
          force(:,iatom) = force(:,iatom) + f_tmp_iatom(:)
       enddo

       !! Add energy coming from sum(KabVaiVbi) term
       e_tmp = 0d0
       do ia_wX= 1,natom_wX_mol(imol)  ! X1,X2,H1,H2
          e_tmp = e_tmp + sum(Kai(ia_wX,:,imol)*Vai(:,imol)) * Vai(ia_wX,imol)
       enddo
       ene_cln = ene_cln - 0.5d0*e_tmp
    enddo      ! imol
!$omp end parallel do

    ! Add force by electric field
    if(flag_use_Et)then
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       force_wX(:,ia_wX,imol) = force_wX(:,ia_wX,imol) + Qai(ia_wX,imol) * Et(:)
    enddo
    enddo
    endif

    ! convert from system with X atoms into ordinary real atom system
    do imol=1,nmol
       !force on H1,H2 (ordinary Coulomb interaction)
       do ia= 2,3  !H1,H2
          iatom = mol2atom_top(imol) + ia-1
          ia_wX = ia+1
          force(:,iatom) = force(:,iatom) + force_wX(:,ia_wX,imol)
       enddo

       ! force through X1,X2 sites in imol & jmol
       do ia= 1,3    ! O,H1,H2 atoms in imol
          iatom = mol2atom_top(imol) + ia-1
          do iX= 1,2       ! X1,X2 sites in imol/jmol
          do ixyz=1,3
             force(ixyz,iatom) = force(ixyz,iatom) + sum(dRXdRatom(:,iX,ixyz,ia,imol)*force_wX(:,iX,imol))
          enddo
          enddo
       enddo
    enddo     ! imol

    ene_vdw = 4d0*(eps_sgm12_vdw_wat * ene_vdw12 - eps_sgm6_vdw_wat * ene_vdw6)
    Uene = Uene + ene_cln + ene_vdw

    !!check of coulomb energy from different calculation
    !e_tmp=0d0
    !do imol=1,nmol
    !do ia_wX=1,natom_wX_mol(imol)
    !   e_tmp = e_tmp + 0.5d0*Qai(ia_wX,imol)*Vai(ia_wX,imol)
    !enddo
    !enddo
    !do imol=1,nmol
    !do ia_wX=1,natom_wX_mol(imol)
    !do ja_wX=1,natom_wX_mol(imol)
    !   e_tmp = e_tmp - 0.5d0*Kai(ia_wX,ja_wX,imol)*Vai(ia_wX,imol)*Vai(ja_wX,imol)
    !enddo
    !enddo
    !enddo
    !write(*,'(a,3e22.10)') "#check-coulomb-energy:", e_tmp, ene_cln, e_tmp-ene_cln

!hoge
!    !write for check 
!    write(*,*) "------ output for check -------"
!    do imol=1,nmol
!    do ia=1,3
!       iatom = mol2atom_top(imol) + ia-1
!       write(*,5) imol,iatom,(Rion(ixyz,iatom),ixyz=1,3)
!    enddo
!    enddo
!  5 format("Rion=",2i4,3e20.10)
!    write(*,10) Uene, (Uene-ene_cln-ene_vdw), ene_vdw, ene_cln
! 10 format("Energy=", 4e20.10)
!    do imol=1,nmol
!       write(*,20) imol,(Q0ai(ia_wX,imol),ia_wX=1,4)
!    enddo
! 20 format("Charge(Q0)=",i4,4e20.10)
!    do imol=1,nmol
!       write(*,25) imol,(Qai(ia_wX,imol),ia_wX=1,4)
!    enddo
! 25 format("Charge(Q)=",i4,4e20.10)
!    do imol=1,nmol
!    do ia=1,3
!       iatom = mol2atom_top(imol) + ia-1
!       write(*,30) imol,iatom,(force(ixyz,iatom),ixyz=1,3)
!    enddo
!    enddo
! 30 format("Force=",2i4,3e20.10)
!    write(*,*) "-------------------------------"
!    stop

  end subroutine force_energy_intermolecular_CRK_cluster

!---------------------------------------------------------

  subroutine image_periodic(drij,idx,nc)
    implicit none
    integer :: j, ipzm,idx(3,nc_max),ipzm_min(3), nc,icx,icy,icz
    real(8) :: tmp_min, drij(3), dist !, drij_new(3)

    if(flag_use_over_cell)then
       !when cut-off radius is larger than half of cell size
       nc = 0
       do icx= -ncell_max(1),ncell_max(1)
       do icy= -ncell_max(2),ncell_max(2)
       do icz= -ncell_max(3),ncell_max(3)
          nc = nc+1
          idx(1,nc) = icx
          idx(2,nc) = icy
          idx(3,nc) = icz
       enddo
       enddo
       enddo

    else
       !when cut-off radius is smaller than half of cell size
       nc=1
       do j=1,3
          tmp_min = 1d99
          do ipzm=-1,1
             dist = abs( drij(j) + ipzm * aL(j) )
             if( dist .le. tmp_min) then
                tmp_min     = dist
                ipzm_min(j) = ipzm
             endif
          enddo
       enddo
       idx(:,nc) = ipzm_min(:)

    endif

  end subroutine

  subroutine get_Xatom_wat(Rxyz_wX,dRXdRatom)
    implicit none
    integer :: i,j, iO,iH1,iH2,iX1,iX2, iatom_O,iatom_H1,iatom_H2, imol
    real(8) :: Rxyz_wX(3,4,nmol), dRXdRatom(3,2,3,3,nmol)
    real(8) :: drOH1_x_drOH2(3,nmol), drOH1(3,nmol), drOH2(3,nmol), rOH1(nmol), rOH2(nmol)
    real(8) :: e_OH1_x_OH2(3), rOH1rOH2, rOH1rOH2_inv, r2_OH1,r2_OH2
    real(8) :: drOH1_x_drOH2_inv, drOH1_x_drOH2_inv2

    iO =1 !(O atom)
    iH1=2 !(H1 atom)
    iH2=3 !(H2 atom)

    iX1=1 !(X1 site)
    iX2=2 !(X2 site)

    do imol=1,nmol

       iatom_O  = mol2atom_top(imol)
       iatom_H1 = mol2atom_top(imol) + 1
       iatom_H2 = mol2atom_top(imol) + 2

       drOH1(:,imol) = Rion(:,iatom_H1) - Rion(:,iatom_O)
       drOH2(:,imol) = Rion(:,iatom_H2) - Rion(:,iatom_O)

       r2_OH1     = sum(drOH1(:,imol)*drOH1(:,imol))
       r2_OH2     = sum(drOH2(:,imol)*drOH2(:,imol))
       rOH1(imol) = sqrt(r2_OH1)
       rOH2(imol) = sqrt(r2_OH2)
       rOH1rOH2     = rOH1(imol) * rOH2(imol)
       rOH1rOH2_inv = 1d0/rOH1rOH2

       drOH1_x_drOH2(1,imol) = drOH1(2,imol)*drOH2(3,imol) - drOH1(3,imol)*drOH2(2,imol)
       drOH1_x_drOH2(2,imol) = drOH1(3,imol)*drOH2(1,imol) - drOH1(1,imol)*drOH2(3,imol)
       drOH1_x_drOH2(3,imol) = drOH1(1,imol)*drOH2(2,imol) - drOH1(2,imol)*drOH2(1,imol)
       drOH1_x_drOH2_inv = 1d0/sqrt(sum(drOH1_x_drOH2(:,imol)*drOH1_x_drOH2(:,imol)))
       drOH1_x_drOH2_inv2= drOH1_x_drOH2_inv * drOH1_x_drOH2_inv
       e_OH1_x_OH2(:) = drOH1_x_drOH2(:,imol) * drOH1_x_drOH2_inv

       Rxyz_wX(:,1,imol) = Rion(:,iatom_O) + dist_OX_wat * e_OH1_x_OH2(:)  ! X1 site
       Rxyz_wX(:,2,imol) = Rion(:,iatom_O) - dist_OX_wat * e_OH1_x_OH2(:)  ! X2 site
       Rxyz_wX(:,3,imol) = Rion(:,iatom_H1)  ! H1 atom
       Rxyz_wX(:,4,imol) = Rion(:,iatom_H2)  ! H2 atom

       ! dR(X1)/dRj(O)
       j=1 !(j=x)
       dRXdRatom(1,iX1,j,iO,imol) = 0d0
       dRXdRatom(2,iX1,j,iO,imol) = ( drOH2(3,imol) - drOH1(3,imol) ) * drOH1_x_drOH2_inv
       dRXdRatom(3,iX1,j,iO,imol) = ( drOH1(2,imol) - drOH2(2,imol) ) * drOH1_x_drOH2_inv
       dRXdRatom(:,iX1,j,iO,imol) = dRXdRatom(:,iX1,j,iO,imol) + e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(2,imol)*(drOH1(3,imol)-drOH2(3,imol)) + drOH1_x_drOH2(3,imol)*(drOH2(2,imol)-drOH1(2,imol)))
       j=2 !(j=y)
       dRXdRatom(1,iX1,j,iO,imol) = ( drOH1(3,imol) - drOH2(3,imol) ) * drOH1_x_drOH2_inv
       dRXdRatom(2,iX1,j,iO,imol) = 0d0
       dRXdRatom(3,iX1,j,iO,imol) = ( drOH2(1,imol) - drOH1(1,imol) ) * drOH1_x_drOH2_inv
       dRXdRatom(:,iX1,j,iO,imol) = dRXdRatom(:,iX1,j,iO,imol) + e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(3,imol)*(drOH1(1,imol)-drOH2(1,imol)) + drOH1_x_drOH2(1,imol)*(drOH2(3,imol)-drOH1(3,imol)))
       j=3 !(j=z)
       dRXdRatom(1,iX1,j,iO,imol) = ( drOH2(2,imol) - drOH1(2,imol) ) * drOH1_x_drOH2_inv
       dRXdRatom(2,iX1,j,iO,imol) = ( drOH1(1,imol) - drOH2(1,imol) ) * drOH1_x_drOH2_inv
       dRXdRatom(3,iX1,j,iO,imol) = 0d0
       dRXdRatom(:,iX1,j,iO,imol) = dRXdRatom(:,iX1,j,iO,imol) + e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(1,imol)*(drOH1(2,imol)-drOH2(2,imol)) + drOH1_x_drOH2(2,imol)*(drOH2(1,imol)-drOH1(1,imol)))

       ! dR(X1)/dRj(H1)
       j=1 !(j=x)
       dRXdRatom(1,iX1,j,iH1,imol) =  0d0
       dRXdRatom(2,iX1,j,iH1,imol) = -drOH2(3,imol) * drOH1_x_drOH2_inv
       dRXdRatom(3,iX1,j,iH1,imol) =  drOH2(2,imol) * drOH1_x_drOH2_inv
       dRXdRatom(:,iX1,j,iH1,imol) = dRXdRatom(:,iX1,j,iH1,imol) + e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(2,imol)*drOH2(3,imol) - drOH1_x_drOH2(3,imol)*drOH2(2,imol))
       j=2 !(j=y)
       dRXdRatom(1,iX1,j,iH1,imol) =  drOH2(3,imol) * drOH1_x_drOH2_inv
       dRXdRatom(2,iX1,j,iH1,imol) =  0d0
       dRXdRatom(3,iX1,j,iH1,imol) = -drOH2(1,imol) * drOH1_x_drOH2_inv
       dRXdRatom(:,iX1,j,iH1,imol) = dRXdRatom(:,iX1,j,iH1,imol) + e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(3,imol)*drOH2(1,imol) - drOH1_x_drOH2(1,imol)*drOH2(3,imol))
       j=3 !(j=z)
       dRXdRatom(1,iX1,j,iH1,imol) = -drOH2(2,imol) * drOH1_x_drOH2_inv
       dRXdRatom(2,iX1,j,iH1,imol) =  drOH2(1,imol) * drOH1_x_drOH2_inv
       dRXdRatom(3,iX1,j,iH1,imol) =  0d0
       dRXdRatom(:,iX1,j,iH1,imol) = dRXdRatom(:,iX1,j,iH1,imol) + e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(1,imol)*drOH2(2,imol) - drOH1_x_drOH2(2,imol)*drOH2(1,imol))

       ! dR(X1)/dRj(H2)
       j=1 !(j=x)
       dRXdRatom(1,iX1,j,iH2,imol) =  0d0
       dRXdRatom(2,iX1,j,iH2,imol) =  drOH1(3,imol) * drOH1_x_drOH2_inv
       dRXdRatom(3,iX1,j,iH2,imol) = -drOH1(2,imol) * drOH1_x_drOH2_inv
       dRXdRatom(:,iX1,j,iH2,imol) = dRXdRatom(:,iX1,j,iH2,imol) - e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(2,imol)*drOH1(3,imol) - drOH1_x_drOH2(3,imol)*drOH1(2,imol))
       j=2 !(j=y)
       dRXdRatom(1,iX1,j,iH2,imol) = -drOH1(3,imol) * drOH1_x_drOH2_inv
       dRXdRatom(2,iX1,j,iH2,imol) =  0d0
       dRXdRatom(3,iX1,j,iH2,imol) =  drOH1(1,imol) * drOH1_x_drOH2_inv
       dRXdRatom(:,iX1,j,iH2,imol) = dRXdRatom(:,iX1,j,iH2,imol) - e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(3,imol)*drOH1(1,imol) - drOH1_x_drOH2(1,imol)*drOH1(3,imol))
       j=3 !(j=z)
       dRXdRatom(1,iX1,j,iH2,imol) =  drOH1(2,imol) * drOH1_x_drOH2_inv
       dRXdRatom(2,iX1,j,iH2,imol) = -drOH1(1,imol) * drOH1_x_drOH2_inv
       dRXdRatom(3,iX1,j,iH2,imol) =  0d0
       dRXdRatom(:,iX1,j,iH2,imol) = dRXdRatom(:,iX1,j,iH2,imol) - e_OH1_x_OH2(:)*drOH1_x_drOH2_inv2 &
         * (drOH1_x_drOH2(1,imol)*drOH1(2,imol) - drOH1_x_drOH2(2,imol)*drOH1(1,imol))

       do i=1,3
       do j=1,3
          dRXdRatom(:,iX1,j,i,imol) =  dRXdRatom(:,iX1,j,i,imol) * dist_OX_wat
          dRXdRatom(:,iX2,j,i,imol) = -dRXdRatom(:,iX1,j,i,imol)
       enddo
       enddo

       do j=1,3
          dRXdRatom(j,iX1,j,iO,imol) = 1d0 + dRXdRatom(j,iX1,j,iO,imol)
          dRXdRatom(j,iX2,j,iO,imol) = 1d0 + dRXdRatom(j,iX2,j,iO,imol)
       enddo

    enddo

  end subroutine


  subroutine get_Q0_K_wat( Q0ai,dQdr_wat,Kai,dKdr_wat,is_wat )
    implicit none
    integer :: imol,ia,iatom,ia_wX,ja_wX,ixyz,is_wat
    real(8) :: Rxyz_imol(3,3), Q0ai(4,nmol),dQ_wat(4,nmol),dQdr_wat(3,4,3,nmol)
    real(8) :: Kai(4,4,nmol),dK_wat(4,4,nmol),dKdr_wat(3,4,4,3,nmol)
    real(8) :: S123(3),dSdr(3,3,3)

    do imol=1,nmol
       do ia=1,natom_mol(imol)
          iatom = mol2atom_top(imol)-1  + ia
          Rxyz_imol(:,ia) = Rion(:,iatom)
       enddo
       call get_conf_dependent_matrix_wat( Rxyz_imol,is_wat,S123,dSdr )

       !correction of Kai and dK/dr = dKdS * dS/dr
       do ia_wX=1,4   !X1,X2,H1,H2
       do ja_wX=1,4   !X1,X2,H1,H2
          dK_wat(ia_wX,ja_wX,imol) = sum(dKdS_wat(ia_wX,ja_wX,:)*S123(:))
          do ia=1,3   !O,H1,H2
          do ixyz=1,3 !x,y,z
             dKdr_wat(ixyz,ia_wX,ja_wX,ia,imol) = sum( dKdS_wat(ia_wX,ja_wX,:)*dSdr(ixyz,:,ia) )
          enddo
          enddo
       enddo
       enddo
       Kai(:,:,imol) = Keq_wat(:,:) + dK_wat(:,:,imol)

       !correction of Q0 and dQ/dr = dQdS * dS/dr
       do ia_wX=1,4   !X1,X2,H1,H2
          dQ_wat(ia_wX,imol) = sum(dQdS_wat(ia_wX,:)*S123(:))
          do ia=1,3   !O,H1,H2
          do ixyz=1,3 !x,y,z
             dQdr_wat(ixyz,ia_wX,ia,imol) = sum( dQdS_wat(ia_wX,:)*dSdr(ixyz,:,ia) )
          enddo
          enddo
       enddo
       Q0ai(:,imol) = Qeq_wat(:) + dQ_wat(:,imol)
    enddo

  end subroutine
  
  subroutine get_conf_dependent_matrix_wat(Rxyz_wat,is_wat,S123,dSdr)
    implicit none
    integer :: is_wat,iatom_O, iatom_H1, iatom_H2, ia
    real(8) :: Rxyz_wat(3,3), rOH1(3),rOH2(3), rOH1_dist,rOH2_dist, drOH1_dist,drOH2_dist
    real(8) :: the_HOH, rOH1_dot_rOH2, cos_the, S123(3),dSdr(3,3,3)
    real(8) :: A, dAdr(3,3), rOH1rOH2, rOH1_dist2, rOH2_dist2

    !calc dr for all 1-2 and 1-3 bonds
    iatom_O  =1
    iatom_H1 =2
    iatom_H2 =3

    rOH1(:) = Rxyz_wat(:,iatom_H1) - Rxyz_wat(:,iatom_O)
    rOH2(:) = Rxyz_wat(:,iatom_H2) - Rxyz_wat(:,iatom_O)
    rOH1_dist  = sqrt(sum(rOH1(:)*rOH1(:)))
    rOH2_dist  = sqrt(sum(rOH2(:)*rOH2(:)))
    drOH1_dist = rOH1_dist - r0_bnd_r(is_wat,1)
    drOH2_dist = rOH2_dist - r0_bnd_r(is_wat,2)

    rOH1_dot_rOH2 = sum(rOH1(:)*rOH2(:))
    cos_the = rOH1_dot_rOH2 / (rOH1_dist * rOH2_dist)
    the_HOH = acos(cos_the)
    
    !internal coordinate
    S123(1) = ( drOH1_dist + drOH2_dist ) / sqrt(2d0)   !=S1(t=1)
    S123(2) = the_HOH - the0_HOH                        !=S2(t=2)
    S123(3) = ( drOH1_dist - drOH2_dist ) / sqrt(2d0)   !=S3(t=3)

    !calculation dSt/dri
    !(for S1 and S3)
    dSdr(:,1,2) =  rOH1(:)/rOH1_dist/sqrt(2d0)   ! dS1/drH1
    dSdr(:,3,2) =  dSdr(:,1,2)                    ! dS3/drH1
    dSdr(:,1,3) =  rOH2(:)/rOH2_dist/sqrt(2d0)   ! dS1/drH2
    dSdr(:,3,3) = -dSdr(:,1,3)                    ! dS3/drH2
    dSdr(:,1,1) = -dSdr(:,1,2) - dSdr(:,1,3)      ! dS1/drO
    dSdr(:,3,1) = -dSdr(:,3,2) - dSdr(:,3,3)      ! dS3/drO
    !(for S2)
    rOH1rOH2   = rOH1_dist * rOH2_dist
    rOH1_dist2 = rOH1_dist * rOH1_dist
    rOH2_dist2 = rOH2_dist * rOH2_dist
    dAdr(:,2) = ( rOH2(:) - rOH1_dot_rOH2 * rOH1(:)/rOH1_dist2 )/rOH1rOH2    !H1
    dAdr(:,3) = ( rOH1(:) - rOH1_dot_rOH2 * rOH2(:)/rOH2_dist2 )/rOH1rOH2    !H2
    dAdr(:,1) = -( rOH1(:)+rOH2(:) - rOH1_dot_rOH2 * ( rOH1(:)/rOH1_dist2 + rOH2(:)/rOH2_dist2 ) )/rOH1rOH2  !O
    A = cos_the
    do ia=1,3   !O,H1,H2
       dSdr(:,2,ia) = -1d0/sqrt(1d0-A*A) * dAdr(:,ia)
    enddo

  end subroutine


  subroutine init_dump_f
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none

    xi_gauss_cln= 0.593d0/ au_length_aa  ![bohr] Gaussian distribution of charges
    if (comm_is_root(nproc_id_global)) &
    write(*,*) "parameter xi for dumping function =", xi_gauss_cln*au_length_aa, "  [A]"

    gm_gauss_cln= 1d0 / sqrt( 2d0*(xi_gauss_cln**2 + xi_gauss_cln**2) )
    coef_derf_gm= 2d0*gm_gauss_cln/sqrt(pi)

    R1_switch   = 4.465d0 / au_length_aa ![bohr] for switching function
    R2_switch   = R1_switch + 1d0        ![bohr]
    R1_switch2  = R1_switch**2
    R2_switch2  = R2_switch**2
    R2_R1_swirtch3_inv = 1d0/( R2_switch2 - R1_switch2 )**3

  end subroutine

  subroutine get_dump_f(dump_f,dump_dfdr,rij,rij2,rij_inv,drij,flag)
    use salmon_math
    implicit none
    integer :: flag
    real(8) :: rij,rij2,rij_inv,drij(3)  !, f_tmp_iatom(3), e_tmp, 
    real(8) :: dump_f, dump_dfdr(3), dump_dfdrij 
    real(8) :: S1_switch,S2_switch, dS1dr_switch,dS2dr_switch, erf_gm_rij
    real(8) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

    if(flag==0)then   !get only dump_f, not the derivative

       if( rij .gt. R2_switch) then  ! no dumping correction for long distance
          dump_f = 1d0
          return
       endif

       erf_gm_rij = erf_salmon(gm_gauss_cln*rij)
       if(     rij.le.R1_switch) then
          dump_f = erf_gm_rij
       else if(rij.le.R2_switch) then
          tmp1 = R2_R1_swirtch3_inv
          tmp2 = ( R2_switch2 - rij2 )**2 * (R2_switch2 + 2d0*rij2 - 3d0*R1_switch2)
          tmp3 = -tmp1
          tmp4 = ( R1_switch2 - rij2 )**2 * (R1_switch2 + 2d0*rij2 - 3d0*R2_switch2)
          S1_switch = tmp2 * tmp1
          S2_switch = tmp4 * tmp3

          tmp5 = R1_switch2 - rij2
          tmp6 = R2_switch2 - rij2
          dS1dr_switch = 12d0*rij * tmp5 * tmp6 * tmp1
          dS2dr_switch = -dS1dr_switch

          dump_f = erf_gm_rij * S1_switch + S2_switch
       else !(for rij.ge.R2_switch)
          dump_f = 1d0
          return
       endif

    else   ! flag \=0, get dump_f and the derivative

       if( rij .gt. R2_switch) then  ! no dumping correction for long distance
          dump_f       = 1d0
          dump_dfdr(:) = 0d0
          return
       endif

       erf_gm_rij = erf_salmon(gm_gauss_cln*rij)
       if(     rij.le.R1_switch) then
          S1_switch   = 1d0
          S2_switch   = 0d0
          dump_dfdrij = coef_derf_gm*exp(-(gm_gauss_cln*rij)**2)
       else if(rij.le.R2_switch) then
          tmp1 = R2_R1_swirtch3_inv
          tmp2 = ( R2_switch2 - rij2 )**2 * (R2_switch2 + 2d0*rij2 - 3d0*R1_switch2)
          tmp3 = -tmp1
          tmp4 = ( R1_switch2 - rij2 )**2 * (R1_switch2 + 2d0*rij2 - 3d0*R2_switch2)
          S1_switch = tmp2 * tmp1
          S2_switch = tmp4 * tmp3
          
          tmp5 = R1_switch2 - rij2
          tmp6 = R2_switch2 - rij2
          dS1dr_switch = 12d0*rij * tmp5 * tmp6 * tmp1
          dS2dr_switch = -dS1dr_switch
          dump_dfdrij  = coef_derf_gm*exp(-(gm_gauss_cln*rij)**2) * S1_switch + erf_gm_rij*dS1dr_switch + dS2dr_switch
       else  !(for rij.ge.R2_switch)
          dump_f       = 1d0
          dump_dfdr(:) = 0d0
          return
       endif
       dump_f       = erf_gm_rij * S1_switch + S2_switch
       dump_dfdr(:) = -dump_dfdrij * rij_inv*drij(:)

    endif !flag

  end subroutine get_dump_f


  subroutine init_ewald_CRK
   !use salmon_global, only: aEwald
    use inputoutput, only: au_length_aa
    implicit none
    integer :: ihx,ihy,ihz, iG, nhx,nhy,nhz
    real(8) :: Gcoef(3), G2_tmp

   !currently alpha of Ewald is from cmd_sc.inp file
   !a_ewald = sqrt(aEwald)  !aEwald from salmon input means (alpha)**2 (not used)

    if (comm_is_root(nproc_id_global)) &
    write(*,*) "alpha of Ewald =", a_ewald/au_length_aa, "  [1/A]"
    coef_derfc = 2d0*a_ewald/sqrt(pi)
    Vuc = aL(1)*aL(2)*aL(3)
    Gcoef(:) = 2d0*pi/aL(:)

    ihx=0
    do
       ihx=ihx+1
       if(abs(Gcoef(1)*ihx).gt.G_cutoff)then
          nhx=ihx
          exit
       endif
    enddo
    ihy=0
    do
       ihy=ihy+1
       if(abs(Gcoef(2)*ihy).gt.G_cutoff)then
          nhy=ihy
          exit
       endif
    enddo
    ihz=0
    do
       ihz=ihz+1
       if(abs(Gcoef(3)*ihz).gt.G_cutoff)then
          nhz=ihz
          exit
       endif
    enddo
    if (comm_is_root(nproc_id_global)) write(*,*) "nhx,nhy,nhz(max)=",nhx,nhy,nhz

    nG_cmd=0
    do ihx= -nhx, nhx
    do ihy= -nhy, nhy
    do ihz= -nhz, nhz
       if(ihx*ihx+ihy*ihy+ihz*ihz==0) cycle  !exclude G=0
       G2_tmp = (Gcoef(1)*ihx)**2 + (Gcoef(2)*ihy)**2 + (Gcoef(3)*ihz)**2
       if(G2_tmp .gt. G_cutoff**2) cycle
       nG_cmd = nG_cmd + 1
    enddo
    enddo
    enddo
    if (comm_is_root(nproc_id_global)) write(*,*) "number of G point in ewald =",nG_cmd
    allocate( Gvec(3,nG_cmd),G2(nG_cmd) )

    iG=0
    do ihx= -nhx, nhx
    do ihy= -nhy, nhy
    do ihz= -nhz, nhz
       if(ihx*ihx+ihy*ihy+ihz*ihz==0) cycle  !exclude G=0
       G2_tmp = (Gcoef(1)*ihx)**2 + (Gcoef(2)*ihy)**2 + (Gcoef(3)*ihz)**2
       if(G2_tmp .gt. G_cutoff**2) cycle
       iG = iG + 1
       Gvec(1,iG) = Gcoef(1)*ihx
       Gvec(2,iG) = Gcoef(2)*ihy
       Gvec(3,iG) = Gcoef(3)*ihz
       G2(iG) = sum(Gvec(:,iG)**2)
    enddo
    enddo
    enddo

    allocate(sumQicosGri(nG_cmd),sumQisinGri(nG_cmd))

  end subroutine

  subroutine get_sumQsinGr_sumQcosGr_ewald(Rxyz_wX,Qai)
    implicit none
    integer :: ia_wX,iG,imol
    real(8) :: GR,Rxyz_wX(3,4,nmol),Qai(4,nmol)

    sumQicosGri(:) = 0d0
    sumQisinGri(:) = 0d0
!$omp parallel do collapse(2) &
!$omp    private(iG,imol,ia_wX,GR) &
!$omp    reduction(+:sumQicosGri,sumQisinGri)
    do iG=1,nG_cmd
    do imol=1,nmol
    do ia_wX=1,natom_wX_mol(imol)
       GR = sum(Gvec(:,iG)*Rxyz_wX(:,ia_wX,imol))
       sumQicosGri(iG) = sumQicosGri(iG) + Qai(ia_wX,imol)*cos(GR)
       sumQisinGri(iG) = sumQisinGri(iG) + Qai(ia_wX,imol)*sin(GR)
    enddo
    enddo
    enddo
!$omp end parallel do

  end subroutine

  
end module

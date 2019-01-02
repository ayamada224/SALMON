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
  !each molecule
  integer,allocatable :: iatom_center(:)
  integer,allocatable :: nbnd_r(:), iatom1_bnd_r(:,:), iatom2_bnd_r(:,:)
  integer,allocatable :: nbnd_a(:), iatom1_bnd_a(:,:), iatom2_bnd_a(:,:)
  real(8),allocatable :: r0_bnd_r(:,:), r0_bnd_a(:,:)
  real(8),allocatable :: kbnd_r(:,:,:,:), kbnd_a(:,:), kbnd_ra(:,:,:), kbnd_rr(:,:,:)
  real(8),allocatable :: charge_cln(:,:), epsilon_vdw(:,:), sigma_vdw(:,:)
  !system
  integer,allocatable :: nb_r, iatom1_b_r(:), iatom2_b_r(:)
  integer,allocatable :: nb_a, iatom1_b_a(:), iatom2_b_a(:)
  integer,allocatable :: nb_rr, ibond1_b_rr(:), ibond2_b_rr(:)
  integer,allocatable :: nb_ra, ibond_b_ra(:),  iangle_b_ra(:)
  real(8),allocatable :: r0_b_r(:), r0_b_a(:)
  real(8),allocatable :: kb_r(:,:,:), kb_a(:), kb_ra(:), kb_rr(:)
  real(8),allocatable :: chg_cln(:), eps_vdw(:), sgm_vdw(:)

  real(8) :: r_cutoff


contains

  subroutine allocate_ff_CRK
    implicit none
    integer :: max_nbnd_r, max_nbnd_a, max_natom

    !change the maximum value if you add larger molecule
    max_natom = 5
    max_nbnd_r  = 5
    max_nbnd_a  = 5

    allocate( natom_mol_s(nmol_s) )    ! # of atoms of the molecule species

    allocate( iatom_center(nmol_s) )
    allocate( nbnd_r(nmol_s), nbnd_a(nmol_s) )
    allocate( iatom1_bnd_r(nmol_s,max_nbnd_r), iatom2_bnd_r(nmol_s,max_nbnd_r) )
    allocate( iatom1_bnd_a(nmol_s,max_nbnd_a), iatom2_bnd_a(nmol_s,max_nbnd_a) )
    allocate( r0_bnd_r(nmol_s,max_nbnd_r), r0_bnd_a(nmol_s,max_nbnd_a) )
    allocate( kbnd_r(nmol_s,max_nbnd_r,2:6,2) )
    allocate( kbnd_a(nmol_s,max_nbnd_a) )
    allocate( kbnd_ra(nmol_s,max_nbnd_r,max_nbnd_a) )
    allocate( kbnd_rr(nmol_s,max_nbnd_r,max_nbnd_r) )

    allocate( charge_cln(nmol_s,max_natom) )
    allocate( epsilon_vdw(nmol_s,max_natom), sigma_vdw(nmol_s,max_natom) )

    iatom1_bnd_r(:,:)=0 ; iatom2_bnd_r(:,:)=0
    iatom1_bnd_a(:,:)=0 ; iatom2_bnd_a(:,:)=0
    r0_bnd_r(:,:)=0d0   ; r0_bnd_a(:,:)=0d0
    kbnd_r(:,:,:,:)=0d0 ; kbnd_a(:,:)=0d0 ; kbnd_ra(:,:,:)=0d0 ; kbnd_rr(:,:,:)=0d0
    charge_cln(:,:)=0d0  ; epsilon_vdw(:,:)=0d0 ; sigma_vdw(:,:)=0d0

  end subroutine

  subroutine load_force_field_CRK(is,name_mol)
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer  is, ib_r,jb_r,ib_a
    real(8) :: the_HOH, Navogadro, eV2J, au2kJmol
    parameter( Navogadro = 6.022140857d23 )
    parameter( eV2J = 1.60218d-19 )
    character(*) :: name_mol

    if(name_mol=="WAT"   .or. &
       name_mol=="WATER" .or. &
       name_mol=="Water" .or. &
       name_mol=="water" ) then
       ! atom order: O, H1, H2

       natom_mol_s(is) = 3
       iatom_center(is) = 1  !Oxigen

       nbnd_r(is) = 2
       nbnd_a(is) = 1

       ! r12 for O-H1 and O-H2 (1-2 interaction)
       do ib_r=1,nbnd_r(is)
          iatom1_bnd_r(is,ib_r) = 1       !O
          iatom2_bnd_r(is,ib_r) = ib_r+1  !H1 or H2

          r0_bnd_r(is,ib_r) = 0.9572d0/au_length_aa ![A]->[au]
          r0_bnd_r(is,ib_r) = 0.9572d0/au_length_aa ![A]->[au]

          kbnd_r(is,ib_r,2,1:2) = (/  0.269d0,  0.269d0 /)  !n=2 [au]
          kbnd_r(is,ib_r,3,1:2) = (/ -0.450d0, -0.500d0 /)  !n=3 [au]
          kbnd_r(is,ib_r,4,1:2) = (/  0.550d0,  3.000d0 /)  !n=4 [au]
          kbnd_r(is,ib_r,5,1:2) = (/ -3.000d0, -3.900d0 /)  !n=5 [au]
          kbnd_r(is,ib_r,6,1:2) = (/  7.000d0, 20.000d0 /)  !n=6 [au]
       enddo

       ! r13 for H-H (1-3 interaction)
       ib_a=1
       iatom1_bnd_a(is,ib_a) = 2  !H1
       iatom2_bnd_a(is,ib_a) = 3  !H2
       the_HOH = 104.52d0 * pi/180d0 ! H-O-H = 104.52 [deg]
       r0_bnd_a(is,ib_a) = 2d0*(0.9572d0/au_length_aa)*sin(the_HOH/2d0)  ![au]
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

       charge_cln(is,1) =-0.329*2  ![e] temporaly 
       charge_cln(is,2) = 0.329
       charge_cln(is,3) = 0.329

       au2kJmol = Navogadro * eV2J / 1d3  * au_energy_ev

       epsilon_vdw(is,1) = 0.6496d0/au2kJmol  ![kJ/mol]->[au]
       epsilon_vdw(is,2) = 0d0
       epsilon_vdw(is,3) = 0d0

       sigma_vdw(is,1) = 3.205d0/au_length_aa  ![A]->[au]
       sigma_vdw(is,2) = 0d0
       sigma_vdw(is,3) = 0d0

    !else if(name_mol=="X") then
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
    integer :: imol, is, iatom, ia
    integer :: ib_r, jb_r, ib_a, iib_r, iib_a, iib_rr, iib_ra, ib_r_last_mol, ib_a_last_mol


    !(calculate nb_r, nb_a)
    nb_r=0
    nb_a=0
    nb_rr=0
    nb_ra=0
    ib_r_last_mol =0
    ib_a_last_mol =0
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
    allocate( iatom1_b_r(nb_r), iatom2_b_r(nb_r) )
    allocate( iatom1_b_a(nb_a), iatom2_b_a(nb_a) )
    allocate( r0_b_r(nb_r), r0_b_a(nb_a) )
    allocate( kb_r(nb_r,2:6,2), kb_a(nb_a) )
    allocate( kb_ra(nb_ra), kb_rr(nb_rr) )
    allocate( ibond1_b_rr(nb_rr), ibond2_b_rr(nb_rr) )
    allocate( ibond_b_ra(nb_ra),  iangle_b_ra(nb_ra) )
    allocate( chg_cln(NI), eps_vdw(NI), sgm_vdw(NI) )


    !(set parameters)
    iib_r  = 0
    iib_a  = 0
    iib_rr = 0
    iib_ra = 0
    do imol=1,nmol
       is = mol2species(imol)
       do ib_r=1,nbnd_r(is)
          iib_r = iib_r + 1
          iatom1_b_r(iib_r) = iatom1_bnd_r(is,ib_r) + mol2atom_top(imol)-1
          iatom2_b_r(iib_r) = iatom2_bnd_r(is,ib_r) + mol2atom_top(imol)-1
          r0_b_r(iib_r) = r0_bnd_r(is,ib_r)
          kb_r(iib_r,2:6,:) = kbnd_r(is,ib_r,2:6,:)
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
          chg_cln(ia) = charge_cln(is,iatom)
          eps_vdw(ia) = epsilon_vdw(is,iatom)
          sgm_vdw(ia) = sigma_vdw(is,iatom)
       enddo

       ib_r_last_mol = ib_r_last_mol + nbnd_r(is)
       ib_a_last_mol = ib_a_last_mol + nbnd_a(is)

    enddo  ! imol

    ! set cutoff radius for non-bonding interaction
    ! currently, half length of the minimum side length of unit cell
    r_cutoff = minval(aL)/2d0
    if (comm_is_root(nproc_id_global)) write(*,*) "cutoff radius =",r_cutoff*au_length_aa," [A]"

  end subroutine


  subroutine cal_force_energy_CRK
    implicit none

    ! Reset forces and total energy
    Uene       = 0d0
    force(:,:) = 0d0

    call force_energy_intramolecular_CRK
    call force_energy_intermolecular_CRK

  end subroutine

  subroutine force_energy_intramolecular_CRK
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer :: ib_r, jb_r, ib_a, ib_rr, ib_ra, iatom1,iatom2,jatom1,jatom2, ipow, isign
    real(8) :: ene1, f_tmp(3), f_tmp_iatom1(3), f_tmp_iatom2(3), f_tmp_jatom1(3), f_tmp_jatom2(3)
    real(8) :: r12_vec(3,nb_r), r12_dist(nb_r), dr12_dist(nb_r), dr12_pow(1:6)
    real(8) :: r13_vec(3,nb_a), r13_dist(nb_a), dr13_dist(nb_a)

    !calc dr for all 1-2 and 1-3 bonds
    do ib_r= 1,nb_r
       iatom1 = iatom1_b_r(ib_r)
       iatom2 = iatom2_b_r(ib_r)
       r12_vec(:,ib_r) = Rion(:,iatom1) - Rion(:,iatom2)
       r12_dist(ib_r) = sqrt(sum(r12_vec(:,ib_r)*r12_vec(:,ib_r)))
       dr12_dist(ib_r) = r12_dist(ib_r) - r0_b_r(ib_r)
    enddo
    do ib_a= 1,nb_a
       iatom1 = iatom1_b_a(ib_a)
       iatom2 = iatom2_b_a(ib_a)
       r13_vec(:,ib_a) = Rion(:,iatom1) - Rion(:,iatom2)
       r13_dist(ib_a) = sqrt(sum(r13_vec(:,ib_a)*r13_vec(:,ib_a)))
       dr13_dist(ib_a) = r13_dist(ib_a) - r0_b_a(ib_a)
    enddo

    ! sum_n kn*(dr12)^n term
    ene1 = 0d0
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

    enddo

    ! krr'*dr12*dr12' term
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

    ! (ktheta,r/2)*(dr13)^2 term
    do ib_a=1,nb_a
       iatom1 = iatom1_b_a(ib_a)
       iatom2 = iatom2_b_a(ib_a)

       f_tmp(:) = -kb_a(ib_a) * dr13_dist(ib_a) * r13_vec(:,ib_a)/r13_dist(ib_a)
       force(:,iatom1) = force(:,iatom1) + f_tmp(:)
       force(:,iatom2) = force(:,iatom2) - f_tmp(:)
       ene1 = ene1 + 0.5d0 * kb_a(ib_a) * dr13_dist(ib_a) * dr13_dist(ib_a)
    enddo

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

    Uene = Uene + ene1

  end subroutine

  subroutine force_energy_intermolecular_CRK
    use inputoutput, only: au_length_aa, au_energy_ev
    implicit none
    integer :: ia,ja, imol, jmol, is, index_image(3)
    integer :: iatom_c, jatom_c, iatom, jatom
    real(8) :: drij(3), rij_inv,rij2,rij2_inv,rij6_inv,rij12_inv,r2_cutoff
    real(8) :: ene_vdw,ene_vdw6,ene_vdw12,ene_cln
    real(8) :: eps_sgm6_vdw_wat, eps_sgm12_vdw_wat
    real(8) :: e_tmp, f_tmp, f_tmp_iatom(3)

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

          ia = mol2atom_top(imol)
          eps_sgm6_vdw_wat = eps_vdw(ia) * sgm_vdw(ia)**6d0
          eps_sgm12_vdw_wat= eps_vdw(ia) * sgm_vdw(ia)**12d0
   
          exit
        endif
    enddo

    r2_cutoff = r_cutoff * r_cutoff


    ene_vdw6  = 0d0
    ene_vdw12 = 0d0
    ene_cln   = 0d0

    do imol=1,nmol
       iatom_c = mol2atom_cnt(imol)
       do jmol=imol+1,nmol
          jatom_c = mol2atom_cnt(jmol)
          drij(:) = Rion(:,jatom_c) - Rion(:,iatom_c)
          call image_periodic(drij,index_image)
          drij(:) = drij(:) + index_image(:)*aL(:)
          rij2 = sum(drij(:)*drij(:))
          if(rij2 .gt. r2_cutoff) cycle   !cut-off

          ! vdw (water, only betwee oxygen atoms)
          rij2_inv = 1d0/rij2
          rij6_inv = rij2_inv*rij2_inv*rij2_inv
          rij12_inv= rij6_inv*rij6_inv
          ene_vdw6  = ene_vdw6  + rij6_inv
          ene_vdw12 = ene_vdw12 + rij12_inv
          f_tmp = 2d0*eps_sgm12_vdw_wat*rij12_inv - eps_sgm6_vdw_wat*rij6_inv
          f_tmp_iatom(:) = -f_tmp * 24d0 * rij2_inv * drij(:)
          force(:,iatom_c) = force(:,iatom_c) + f_tmp_iatom(:)
          force(:,jatom_c) = force(:,jatom_c) - f_tmp_iatom(:)

          ! coulomb (currently, simple coulomb interaction for test)
          do ia= 1,natom_mol(imol)
             iatom = mol2atom_top(imol) + ia-1
             do ja= 1,natom_mol(jmol)
                jatom = mol2atom_top(jmol) + ja-1

                drij(:) = (Rion(:,jatom) - Rion(:,iatom)) + index_image(:)*aL(:)
                rij2    = sum(drij(:)*drij(:))
                rij2_inv= 1d0/rij2
                rij_inv = sqrt(rij2_inv)

                e_tmp = chg_cln(iatom)*chg_cln(jatom) * rij_inv
                f_tmp_iatom(:) = -e_tmp * rij2_inv * drij(:)
                force(:,iatom) = force(:,iatom) + f_tmp_iatom(:)
                force(:,jatom) = force(:,jatom) - f_tmp_iatom(:)
                ene_cln   = ene_cln + e_tmp
             enddo
          enddo

       enddo
    enddo

    ene_vdw = 4d0*(eps_sgm12_vdw_wat * ene_vdw12 - eps_sgm6_vdw_wat * ene_vdw6)
    Uene = Uene + ene_cln + ene_vdw

  end subroutine

  subroutine image_periodic(drij,ipzm_min)
    implicit none
    integer :: j, ipzm, ipzm_min(3)
    real(8) :: tmp_min, drij(3), dist !, drij_new(3)

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

  end subroutine

end module

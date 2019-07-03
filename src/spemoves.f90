module spemoves
!! Additional MC moves (swap etc..).
use state
use moves
implicit none

contains

subroutine init_moves(st, mvs)
  !! Initialize all builtin moves.
  type(PcmcState) :: st
  type(MCMoves) :: mvs
  real(8) :: prob_volmod
  call moves_reset(mvs)
  call add_move(mvs, "translations", 1.d0, translation1P)
  call add_move(mvs,"swaps"                                           , st%inp%prob_swap,echang2p)
  call add_move(mvs,"interbox swaps"                                  , st%inp%prob_swap_box,echang2p_box)
  call add_move(mvs,"Particule exchange"                              , st%inp%prob_box2box,mouv1p_box)
  call add_move(mvs,"volume exchange"                                 , st%inp%prob_volex,echangvol)
  call add_move(mvs,"vol + part exchange"                             , st%inp%prob_volparex,echangvol_o_mouv1p_box)
  call add_move(mvs,"part exchange rosenbluthed"                      , st%inp%prob_b2b_rosbth,mouv1p_box_rosbth)
  call add_move(mvs,"interbox swaps impur"                            , st%inp%prob_swap_impur , echang2p_box_impur)
  call add_move(mvs,"interbox swaps impur rosenbluthed"               , st%inp%prob_swap_impur_rosbth , &
    echang2p_box_impur_rosbth)
  call add_move(mvs,"interbox swaps impur rosenbluthed 2"             , st%inp%prob_swap_impur_rosbth2 , &
    echang2p_box_impur_rosbth2)
  if(st%inp%asimu(1:3)=="NPT") then ! NVT + Gibbs by settings probs is def
      if(st%inp%prob_volmod<=0.d0) then
        prob_volmod=1.d0/dble(maxval(st%boxes(1:st%ntotbox)%n))
      else
        prob_volmod = st%inp%prob_volmod
      endif
      call add_move(mvs,"volume changes"            , prob_volmod , changevol_npt)
  endif
end subroutine init_moves

integer function echang2p(st, nr)
  !! Attempts to swap two random particles in a box  
  use swapmap
  implicit none
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  type(NrgRec) :: ener1,ener0
  logical :: bool,authtoselfam
  integer :: i,ip,ip1,ibox,fam1,counter=0
  counter=counter+1
  ibox = 1 + floor(ran2(st%rng)*st%ntotbox)
  ener1 = st%boxes(ibox)%ener
  ener0=ener1
  
  ip = 1 + floor(ran2(st%rng)*st%boxes(ibox)%n)
  ip1 = 1 + floor(ran2(st%rng)*st%boxes(ibox)%n)
  fam1 = st%boxes(ibox)%parts(ip)%Famille
  authtoselfam = .false.
  do i=1,st%dist%nfam
     if(i/=fam1 .and. st%boxes(ibox)%fam(i)>0) then
       authtoselfam = .true.
       exit
     end if
  end do
  if(authtoselfam) then
     do while( st%boxes(ibox)%parts(ip1)%famille == fam1)
        ip1 = 1 + floor(ran2(st%rng)*st%boxes(ibox)%n)
     enddo
  else
     echang2p = mv_cir_rej
     return
  end if

  associate(labox => st%boxes(ibox),part1 => st%boxes(ibox)%parts(ip),part2 => st%boxes(ibox)%parts(ip1))

    ener1 = ener0 - nr%energy_1p(st,labox,part1) - nr%energy_1p(st,labox,part2)

    call swap_2parts(st,ibox, part1, part2)
    ener1 = ener1 + nr%energy_1p(st,labox,part1)
    if(ener1%noverlap) ener1 = ener1 + nr%energy_1p(st,labox,part2)
    if(ener1%noverlap) then
      if(metropolis((ener1%ep-ener0%ep)/kT,st%rng)) then
          labox%ener=ener1
          echang2p=mv_acc
          call swapmap_update(part1%famille,part2%famille)
       else
          call swap_2parts(st,ibox, part1, part2)
          echang2p=mv_nrj_rej
       endif
    else
       call swap_2parts(st,ibox, part1, part2)
       echang2p=mv_hs_rej
    endif
  end associate
end function echang2p

integer function changevol_npt(st, nr) result(ret)
  !! Volume change in the isobaric ensemble.
  implicit none
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  type(NrgRec) :: ener1_1,ener0_1
  integer :: ibox1
  real(8) :: deltabmu, Pref
  real(8) :: fac1,dvol,rcell1_0(3),dumpv(3),oldbox(3,3)


  Pref = st%inp%Pref / (unitPsT*st%inp%temp) ! convert to kT units

  ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
  
  associate(labox1 => st%boxes(ibox1))

    ener0_1 = labox1%ener

    dvol = (ran2(st%rng)-0.5d0)*st%inp%dv_vol
    if(st%inp%std_gibbs_volex) then
       fac1=1.d0 + dvol/labox1%volume
    else
       fac1 = exp(dvol)
    end if
    if(fac1<=0.d0) then    
       ret=mv_cir_rej
       return
    end if
    if(st%inp%std_gibbs_volex) then
       deltabmu = - (labox1%n)*log(fac1)  
    else
       deltabmu = - (labox1%n+1)*log(fac1)
    end if
    deltabmu = deltabmu + Pref*labox1%volume*(fac1-1.d0)
    oldbox = labox1%tens
    call box_set_tensor(labox1, fac1**(1.d0/3.d0)*labox1%tens)
    if(st%celllist) then
       rcell1_0 = st%clistdecomp(ibox1)%rcell
       dumpv(1)= st%rc/(labox1%tens(1,1))
       dumpv(2)= st%rc/(labox1%tens(2,2))
       dumpv(3)= st%rc/(labox1%tens(3,3))
       call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,&
            dumpv)
    endif
    ener1_1 = nr%energy_box(st,labox1)
    if(ener1_1%noverlap) then
       if(metropolis((ener1_1%ep-ener0_1%ep)/kT + deltabmu,st%rng)) then
          labox1%ener = ener1_1
          ret=mv_acc
       else
          call box_set_tensor(labox1, oldbox)
          if(st%celllist) then
            call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,rcell1_0)
          end if
          ret=mv_nrj_rej
       end if
    else
      call box_set_tensor(labox1, oldbox)
      if(st%celllist) then
        call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,rcell1_0)
      end if
      ret=mv_hs_rej
    end if
  end associate

end function changevol_npt

integer function echang2p_box(st, nr)
  !! attempts to swap two particles located in two different boxes
  use swapmap
  implicit none
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2
  logical :: bool,authtoselfam
  integer :: i,ip,ip1,ibox1,ibox2,fam1,fam2,idump,idump1
  real(8) :: deltabmu
  ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
  ip = 1 + floor(ran2(st%rng)*st%boxes(ibox1)%n)
  fam1 = st%boxes(ibox1)%parts(ip)%famille

  ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
  if(st%ntotbox>1) then
     do while(ibox1 == ibox2)
        ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
     end do
  else
     echang2p_box=mv_cir_rej
     return
  end if
  
  associate(labox1 => st%boxes(ibox1), labox2 => st%boxes(ibox2))
  
  if(labox1%n<1.or.labox2%n<1) then
     echang2p_box=mv_cir_rej
     return
  end if
  ip1 = 1 + floor(ran2(st%rng)*labox2%n)
  authtoselfam = .false.
  do i=1,st%dist%nfam
     if(i/=fam1 .and. labox2%fam(i)>0) authtoselfam = .true.
     exit
  end do
  if(authtoselfam) then
     do while(labox2%Parts(ip1)%Famille == fam1 )
        ip1 = 1 + floor(ran2(st%rng)*labox2%n)
     enddo
  else
     echang2p_box = mv_cir_rej
     return
  end if

  ener1_1 = labox1%ener
  ener1_2 = labox2%ener
  ener0_1=ener1_1
  ener0_2=ener1_2
  
  associate(part1 => labox1%parts(ip), part2 => labox2%parts(ip1))

  fam1=part1%famille
  fam2=part2%famille

  idump=max(labox1%fam(fam1),1)
  idump1=max(labox2%fam(fam2),1)
  deltabmu = log(((1.d0*(labox1%fam(fam2)+1))*(1.d0*(labox2%fam(fam1)+1)))/&
       ((1.d0*idump)*(1.d0*idump1)))
  
  ener1_1 = ener0_1 - nr%energy_1p(st,labox1,part1)
  ener1_2 = ener0_2 - nr%energy_1p(st,labox2,part2)

  call swap_2parts_2box(st,ibox1, ibox2, part1, part2)
  ener1_1 = ener1_1 + nr%energy_1p(st,labox1,part1)
  if(ener1_1%noverlap) ener1_2 = ener1_2 + nr%energy_1p(st,labox2,part2)
  if(ener1_1%noverlap .and. ener1_2%noverlap) then
     if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT + deltabmu,st%rng)) then
        labox1%ener=ener1_1
        labox2%ener=ener1_2
        echang2p_box=mv_acc
        call swapmap_box_update(part1%famille,part2%famille)
     else
        call swap_2parts_2box(st,ibox2, ibox1, part2, part1)
        echang2p_box=mv_nrj_rej
     endif
  else
     call swap_2parts_2box(st,ibox2, ibox1, part2, part1)
     echang2p_box=mv_hs_rej
  endif

  end associate
  end associate
end function echang2p_box

  integer function mouv1p_box(st, nr)
    !! Attempts to move a particle from a box to a different box at a random position.
    use iso_fortran_env
    use chempot
    implicit none
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2
    logical :: bool
    integer :: ip,ibox1,ibox2,i,idump,idump1,fam,counter=0
    real(8) :: store_pos(3)
    real(8) :: deltabmu
    counter=counter+1
    ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
    
    if(st%boxes(ibox1)%n<2) then
       mouv1p_box=mv_cir_rej
       return
    end if
    ip = 1 + floor(ran2(st%rng)*st%boxes(ibox1)%n)
    

    ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
    if(st%ntotbox>1) then
       do while(ibox1 == ibox2)
          ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
       end do
    else
       mouv1p_box=mv_cir_rej
       return
    end if

    associate(part => st%boxes(ibox1)%Parts(ip), labox1 => st%boxes(ibox1),labox2 => st%boxes(ibox2))

    ! check particle vector
    if(labox2%n >= st%Nbuff) then
       write(error_unit,'(A,I0,A,I0,A)') "MOUV1P_BOX :: Particle vector of box ",ibox2&
            ," almost exploded (buffer size = ",st%Nbuff,"). Bye."
       write(*,*)
       stop 3
    end if
    
    ener1_1 = labox1%ener
    ener1_2 = labox2%ener
    ener0_1=ener1_1
    ener0_2=ener1_2

    store_pos=part%pos
    ener1_1 = ener0_1 - nr%energy_1p(st,labox1,part)

    do i=1,3
       part%pos(i) = ran2(st%rng)
    end do
    !
    fam = part%famille
    idump=max(labox1%fam(fam),1)
    idump1=max(labox2%fam(fam),1)
    deltabmu = - log((1.d0*idump)/(1.d0*(idump1+1))*labox2%volume/labox1%volume)

    ener1_2 =  ener0_2 + nr%energy_1p(st,labox2,part)
    if(ener1_2%noverlap) then
       if(chempwid_doit_doit) then
          chempwid_bmu(part%famille,ibox2) = chempwid_bmu(part%famille,ibox2) + &
               exp(ener0_2%ep - ener1_2%ep)
       end if
       if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT + deltabmu,st%rng)) then
          labox1%ener=ener1_1
          labox2%ener=ener1_2
          call transfer_part_b2b(st,part,ibox1,ibox2)
          mouv1p_box=mv_acc
       else
          part%pos=store_pos
          mouv1p_box=mv_nrj_rej
       endif
    else
       part%pos=store_pos
       mouv1p_box=mv_hs_rej
    endif
    if(chempwid_doit_doit) chempwid_nsample(part%famille,ibox2) = chempwid_nsample(part%famille,ibox2) + 1

  end associate
  end function mouv1p_box

  integer function echangvol(st, nr)
    !! Attempts to exchange a quantity of volume between two boxes
    implicit none
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2
    integer :: ibox1,ibox2,counter=0
    real(8) :: deltabmu
    real(8) :: fac1,fac2,dvol,vsum,rcell1_0(3),rcell2_0(3),dumpv(3)
    real(8) :: oldbox1(3,3), oldbox2(3,3)
    counter=counter+1
    ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)

    ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)

    if(st%ntotbox>1) then
       do while(ibox1 == ibox2)
          ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
       end do
    else
       ECHANGVOL=mv_cir_rej
       return
    end if

    associate(labox1 => st%boxes(ibox1), labox2=> st%boxes(ibox2))

    ener0_1 = labox1%ener
    ener0_2 = labox2%ener


    dvol = ran2(st%rng)*st%inp%dv_volex
    if(st%inp%std_gibbs_volex) then
       fac1=1.D0 + dvol/labox1%volume
       fac2=1.D0 - dvol/labox2%volume
    else
       dvol = exp(dvol)
       vsum = labox1%volume+labox2%volume
       fac1=(dvol*labox1%volume*vsum)/(vsum+labox1%volume*(dvol-1.D0))  !Vnew1
       fac2=vsum-fac1  ! Vnew2
       fac1=fac1/labox1%volume
       fac2=fac2/labox2%volume
    end if
    if(fac1<=0.D0 .or. fac2<=0.D0) then    
       echangvol=mv_cir_rej
       return
    end if
    if(st%inp%std_gibbs_volex) then
       deltabmu = - (labox1%n)*log(fac1) - (labox2%n)*log(fac2)
    else
       deltabmu = - (labox1%n+1)*log(fac1) - (labox2%n+1)*log(fac2)
    end if
    oldbox1 = labox1%tens
    call box_set_tensor(labox1, labox1%tens*(fac1)**(1.D0/3.D0))
    oldbox2 = labox2%tens
    call box_set_tensor(labox2,  labox2%tens*(fac2)**(1.D0/3.D0))
    if(st%celllist) then
       rcell1_0 = st%clistdecomp(ibox1)%rcell
       dumpv(1)= st%rc/(labox1%tens(1,1)*(fac1)**(1.D0/3.D0))
       dumpv(2)= st%rc/(labox1%tens(2,2)*(fac1)**(1.D0/3.D0))
       dumpv(3)= st%rc/(labox1%tens(3,3)*(fac1)**(1.D0/3.D0))
       call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,&
            dumpv)
       rcell2_0 = st%clistdecomp(ibox2)%rcell
       dumpv(1)= st%rc/(labox2%tens(1,1)*(fac2)**(1.D0/3.D0))
       dumpv(2)= st%rc/(labox2%tens(2,2)*(fac2)**(1.D0/3.D0))
       dumpv(3)= st%rc/(labox2%tens(3,3)*(fac2)**(1.D0/3.D0))
       call celldec_cree(st%clistdecomp(ibox2),labox2%n,labox2%parts,&
            dumpv)
    endif
    ener1_1 = nr%energy_box(st,labox1)
    if(ener1_1%noverlap) then
       ener1_2 = nr%energy_box(st,labox2)
       if(ener1_2%noverlap) then
          if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT + deltabmu,st%rng)) then
             labox1%ener = ener1_1
             labox2%ener = ener1_2
             echangvol=mv_acc
          else
            call rejection()
            echangvol=mv_nrj_rej
          end if
       else
          call rejection()
          echangvol=mv_hs_rej
       end if
    else
      call rejection()
      echangvol=mv_hs_rej
    end if

    end associate

    contains

      subroutine rejection()
        associate(labox1 => st%boxes(ibox1), labox2=> st%boxes(ibox2))
        call box_set_tensor(labox1, oldbox1)
        call box_set_tensor(labox2, oldbox2)
        if(st%celllist) then
          call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,rcell1_0)
          call celldec_cree(st%clistdecomp(ibox2),labox2%n,labox2%parts,rcell2_0)
        end if

        end associate
      end subroutine
  
  end function echangvol

  integer function echangvol_o_mouv1p_box(st, nr)
    !! Attempts to exchange a quantity of volume between two boxes,
    !! and move a particule from one box to the other in a single move.
    !!
    !! This was implement hoping to improve the acceptance of dense systems
    !! but it does not seem to have a significant impact. In short, this move is
    !! probably useless.    
    implicit none
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2
    integer :: i,ip,ibox1,ibox2,fam,idump,idump1
    real(8) :: deltabmu
    real(8) :: fac1,fac2,dvol,vsum,rcell1_0(3),rcell2_0(3)
    real(8) :: oldbox1(3,3), oldbox2(3,3)
    type(Particule) :: part
    
    ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
    if(st%boxes(ibox1)%n<2) then
       echangvol_o_mouv1p_box=mv_cir_rej
       return
    end if
    ip = 1 + floor(ran2(st%rng)*st%boxes(ibox1)%n)

    ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
    if(st%ntotbox>1) then
       do while(ibox1 == ibox2)
          ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
       end do
    else
       echangvol_o_mouv1p_box=mv_cir_rej
       return
    end if

    associate(labox1 => st%boxes(ibox1), labox2=> st%boxes(ibox2))

    ener0_1 = labox1%ener
    ener0_2 = labox2%ener

    dvol = (2.d0*ran2(st%rng)-1.d0)*st%inp%dv_volex
    if(st%inp%std_gibbs_volex) then
       fac1=1.d0 + dvol/labox1%volume
       fac2=1.d0 - dvol/labox2%volume
    else
       dvol = exp(dvol)
       vsum = labox1%volume+labox2%volume
       fac1=(dvol*labox1%volume*vsum)/(vsum+labox1%volume*(dvol-1.D0))
       fac2=vsum-fac1
       fac1=fac1/labox1%volume
       fac2=fac2/labox2%volume
    end if
    if(fac1<=0.d0 .or. fac2<=0.d0) then
       echangvol_o_mouv1p_box=mv_cir_rej
       return
    end if
    if(st%inp%std_gibbs_volex) then
       deltabmu = - (labox1%n)*log(fac1) - (labox2%n)*log(fac2)
    else
       deltabmu = - (labox1%n+1)*log(fac1) - (labox2%n+1)*log(fac2)
    endif

    part = labox1%parts(ip)
    call add_part_box_nocelllist(st,ibox2,part)   
    call remove_part_box_nocelllist(st,ibox1,ip, part)
    idump=labox2%n
    do i=1,3
       labox2%parts(idump)%pos(i) = ran2(st%rng)
    end do
    fam = part%famille
    idump=max(labox1%fam(fam),1)
    idump1=max(labox2%fam(fam),1)
    deltabmu = deltabmu - log((1.D0*idump)/(1.D0*(idump1+1))*(fac2*labox2%volume)/(fac1*labox1%volume))
    oldbox1 = labox1%tens
    call box_set_tensor(labox1, labox1%tens*(fac1)**(1.D0/3.D0))
    oldbox2 = labox2%tens
    call box_set_tensor(labox2, labox2%tens*(fac2)**(1.D0/3.D0))
    if(st%celllist) then
       rcell1_0 = st%clistdecomp(ibox1)%rcell
       call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,&
            st%clistdecomp(ibox1)%rcell*(fac1)**(1.D0/3.D0))
       rcell2_0 = st%clistdecomp(ibox2)%rcell
       call celldec_cree(st%clistdecomp(ibox2),labox2%n,labox2%parts,&
            st%clistdecomp(ibox2)%rcell*(fac2)**(1.D0/3.D0))
    endif
    ener1_1 = nr%energy_box(st,labox1)
    if(ener1_1%noverlap) then
       ener1_2 = nr%energy_box(st,labox2)
       if(ener1_2%noverlap) then
          if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT + deltabmu,st%rng)) then
             labox1%ener = ener1_1
             labox2%ener = ener1_2
             echangvol_o_mouv1p_box=mv_acc
          else
            call rejection()
             echangvol_o_mouv1p_box=mv_nrj_rej
          end if
       else
          call rejection()
          echangvol_o_mouv1p_box=mv_hs_rej
       end if
    else
      call rejection()
      echangvol_o_mouv1p_box=mv_hs_rej
    end if

    end associate

    contains

    subroutine rejection()
      associate(labox1 => st%boxes(ibox1), labox2=> st%boxes(ibox2))

      call add_part_box_nocelllist(st,ibox1,part)
      idump=labox2%n
      call remove_part_box_nocelllist(st,ibox2,idump, part)
      call box_set_tensor(labox1, oldbox1)
      call box_set_tensor(labox2, oldbox2)
      if(st%celllist) then
        call celldec_cree(st%clistdecomp(ibox1),labox1%n,labox1%parts,rcell1_0)
        call celldec_cree(st%clistdecomp(ibox2),labox2%n,labox2%parts,rcell2_0)
      end if

      end associate
    end subroutine
        
  end function echangvol_o_mouv1p_box

  

  integer function mouv1p_box_rosbth(st, nr)
    !! Attempts to move a particle to a different box with rosenbluth trials
    !!
    !! A series of trial positions in the new box is generated, each with a weight equal to
    !! \( w= \exp(-\beta E(trial)-E(old)) \). The finally attempted position is selected
    !! at random according to the weights. As this procedure introduce a bias,
    !! weights of  the old configuration and weights of the new are computed
    !! and injected into the acceptance ratio to correct the bias.
    use iso_fortran_env
    use chempot
    implicit none
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2,ener_dump
    logical :: bool,authtoselfam
    integer :: ip,ibox1,ibox2,i,idump,idump1,fam,k,k_selected,counter=0,counter2=0
    integer :: ntrials_rosbth
    real(8) :: store_pos(3),dump
    real(8),allocatable :: trial_pos(:,:),trial_prob(:)
    real(8) :: deltabmu, rosbth_fact,weight_new,weight_old,dumppos(3)
    type(NrgRec), allocatable :: Trial_ener(:)
    save ::  trial_pos,trial_prob, Trial_ener
    counter=counter+1
    ntrials_rosbth = st%inp%ntrials_rosbth
    if(.not.allocated(trial_pos)) then
       allocate(trial_pos(3,ntrials_rosbth),trial_prob(ntrials_rosbth),Trial_ener(ntrials_rosbth))
       counter2=counter2+1
    end if
    ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
    associate(labox1 => st%boxes(ibox1))
    if(labox1%n<2) then
       mouv1p_box_rosbth=mv_cir_rej
       return
    end if
    ip = 1 + floor(ran2(st%rng)*labox1%n)
    authtoselfam = .false.
    do i=1,st%inp%n_rospechmax
       if(labox1%fam(i)>0) authtoselfam = .true.
       exit
    end do
    if(authtoselfam) then
       do while(labox1%parts(ip)%famille > st%inp%n_rospechmax)
          ip = 1 + floor(ran2(st%rng)*labox1%n)
       end do
    end if

    ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
    if(st%ntotbox>1) then
       do while(ibox1 == ibox2)
          ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
       end do
    else
       mouv1p_box_rosbth=mv_cir_rej
       return
    end if

    associate(part => labox1%Parts(ip), labox2 => st%boxes(ibox2))

    ! check particle vector
    if(labox2%n >= st%Nbuff) then
       write(error_unit,'(A,I0,A,I0,A)') "MOUV1P_BOX :: Particle vector of box ",ibox2&
            ," almost exploded (buffer size = ",st%Nbuff,"). Bye."
       write(*,*)
       stop 3
    end if
    
    ener1_1 = labox1%ener
    ener1_2 = labox2%ener
    ener0_1=ener1_1
    ener0_2=ener1_2

    store_pos=part%pos
    ener_dump = nr%energy_1p(st,labox1,part)
    ener1_1 = ener0_1 - ener_dump

    ! Old trials
    weight_old = exp(-ener_dump%ep/kT)
    do k=2,ntrials_rosbth
       do i=1,3
          dumppos(i) = ran2(st%rng)
       end do
       part%pos(:) = dumppos(:)
       ener_dump = nr%energy_1p(st,labox1,part)
       if(ener_dump%noverlap) then
          dump = exp(-ener_dump%ep/kT)
       else
          dump = 0.D0
       endif
       weight_old = weight_old +  dump
    end do

    ! "New" trials
    weight_new = 0.D0
    do k=1,ntrials_rosbth
       do i=1,3
          trial_pos(i,k) = ran2(st%rng)
       end do
       part%pos(:) = trial_pos(1:3,k)
       ener_dump = nr%energy_1p(st,labox2,part)
       if(ener_dump%noverlap) then
          trial_prob(k) = exp(-ener_dump%ep/kT)
          trial_ener(k) = ener_dump
       else
          trial_prob(k) = 0.D0
       endif
       weight_new = weight_new +  trial_prob(k)
    end do
    if(chempwid_doit_doit) then
       chempros_bmu(part%famille,ibox2) = chempros_bmu(part%famille,ibox2) + weight_new/ntrials_rosbth
       chempros_nsample(part%famille,ibox2) = chempros_nsample(part%famille,ibox2) + 1
    end if
    if(weight_new < 1.D-13) then ! Arbitrary thing Warning
       part%pos = store_pos
       if(weight_new==0.D0) then
          mouv1p_box_rosbth=mv_hs_rej
       else
          mouv1p_box_rosbth=mv_nrj_rej
       end if
       return
    end if
    do k=2,ntrials_rosbth
       trial_prob(k)= trial_prob(k)+trial_prob(k-1)
    end do
    do k=1,ntrials_rosbth
       trial_prob(k) = trial_prob(k)/weight_new
    end do
    dump = ran2(st%rng)
    k_selected = ntrials_rosbth
    do k=1,ntrials_rosbth
       if(dump<=trial_prob(k)) then
          k_selected = k
          exit
       end if
    end do
    if(k_selected < 1 .or. k_selected > ntrials_rosbth) then
       write(*,*) real(trial_prob)
       stop 666
    end if
    part%pos = trial_pos(1:3,k_selected)
    ener1_2%ep = ener0_2%ep + Trial_ener(k_selected)%ep
    rosbth_fact = weight_new/weight_old
    
    fam = part%famille
    idump=max(labox1%fam(fam),1)
    idump1=max(labox2%fam(fam),1)
    deltabmu = - log((1.D0*idump)/(1.D0*(idump1+1))*labox2%volume/labox1%volume)

    if(trial_ener(k_selected)%noverlap) then
       if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT  + deltabmu - log(rosbth_fact),st%rng)) then
          labox1%ener=ener1_1
          labox2%ener=ener1_2
          call transfer_part_b2b(st,part,ibox1,ibox2)
          mouv1p_box_rosbth=mv_acc
       else
          part%pos=store_pos
          mouv1p_box_rosbth=mv_nrj_rej
       endif
    else
       part%pos=store_pos
       mouv1p_box_rosbth=mv_hs_rej
    endif

    end associate
    end associate

  end function mouv1p_box_rosbth

  

  integer function echang2p_box_impur(st, nr)
    !! Attempts to swap a particle with an "impurity"
    !! (ie the first particle family of particle, supposedly very small)
    !! between two different boxes.
    !!
    !! This was introduced to accelerate the equilibration of dense systems.
    !! Indeed, a very small particle (in small quantities, this can be taken as
    !! as an impurity) can move more easily from one box to another.
    !! Once there, it can help bigger particles to pass to a different box by swapping
    !! their positions. In pratice, it did not help in our case.
    use swapmap
    implicit none
    type(PcmcState):: st
    type(NrgRoutines) :: nr
    type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2
    logical :: bool,authtoselfam
    integer :: ip,ip1,ibox1,ibox2,fam1,fam2,idump,idump1
    real(8) :: deltabmu
    type(Particule) :: store_part
  
    ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
    associate(labox1 => st%boxes(ibox1))
    ip = 1 + floor(ran2(st%rng)*labox1%n)
    fam1 = labox1%parts(ip)%famille
    if(all(labox1%fam(2:st%dist%nfam)==0)) then
      echang2p_box_impur=mv_cir_rej
      return
    end if
    do while(fam1==1)
      ip = 1 + floor(ran2(st%rng)*labox1%n)
      fam1 = labox1%parts(ip)%famille
    end do
  
    ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
    if(st%ntotbox>1) then
      do while(ibox1 == ibox2)
        ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
      end do
    else
      echang2p_box_impur=mv_cir_rej
      return
    end if
  
    associate(labox2 => st%boxes(ibox2))
    if(labox1%n<1.or.labox2%n<1) then
      echang2p_box_impur=mv_cir_rej
      return
    end if
    ip1 = 1 + floor(ran2(st%rng)*labox2%n)

    if(fam1/=1 .and. labox2%fam(1)>0) authtoselfam = .true.        ! impurs are supposed to be fam 1
    if(authtoselfam) then
     do while(labox2%Parts(ip1)%Famille /= 1 )
       ip1 = 1 + floor(ran2(st%rng)*labox2%n)
     enddo
   else
     echang2p_box_impur=mv_cir_rej
     return
   end if
  
   ener1_1 = labox1%ener
   ener1_2 = labox2%ener
   ener0_1=ener1_1
   ener0_2=ener1_2
  
   associate(part1 => labox1%parts(ip), part2 => labox2%parts(ip1))
   fam1=part1%famille
   fam2=part2%famille

   idump=max(labox1%fam(fam1),1)
   idump1=max(labox2%fam(fam2),1)
   deltabmu = log(((1.d0*(labox1%fam(fam2)+1))*(1.d0*(labox2%fam(fam1)+1)))/&
        ((1.d0*idump)*(1.d0*idump1)))
  
   ener1_1 = ener1_1 - nr%energy_1p(st,labox1,part1)
   ener1_2 = ener1_2 - nr%energy_1p(st,labox2,part2)

   store_part = part1
   part1%famille=part2%famille
   part1%numero=part2%numero
   part1%rayon=part2%rayon
   part2%famille=store_part%famille
   part2%numero = store_part%numero
   part2%rayon = store_part%rayon
   ener1_1 = ener1_1 + nr%energy_1p(st,labox1,part1)
   if(ener1_1%noverlap) ener1_2 = ener1_2 +  nr%energy_1p(st,labox2,part2)
   if(ener1_1%noverlap .and. ener1_2%noverlap) then
     if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT + deltabmu,st%rng)) then
       labox1%ener=ener1_1
       labox2%ener=ener1_2
       labox1%fam(fam1)=labox1%fam(fam1)-1
       labox1%fam(fam2)=labox1%fam(fam2)+1
       labox2%fam(fam2)=labox2%fam(fam2)-1
       labox2%fam(fam1)=labox2%fam(fam1)+1
       echang2p_box_impur=mv_acc
       call swapmap_impur_update(store_part%famille)
     else
       part2%famille=part1%famille
       part2%numero=part1%numero
       part2%rayon=part1%rayon
       part1%famille=store_part%famille
       part1%numero = store_part%numero
       part1%rayon = store_part%rayon
       echang2p_box_impur=mv_nrj_rej
     endif
   else
     part2%famille=part1%famille
     part2%numero=part1%numero
     part2%rayon=part1%rayon
     part1%famille=store_part%famille
     part1%numero = store_part%numero
     part1%rayon = store_part%rayon
     echang2p_box_impur=mv_hs_rej
   endif

   end associate
   end associate
   end associate
 end function echang2p_box_impur

 integer function echang2p_box_impur_rosbth(st, nr)
  !! Attempts to swap a particle with an "impurity"
  !! (ie the first particle family of particle, supposedly very small)
  !! between two different boxes, with rosenbluth trials.
  !!
  !! Same as \ref echang2p_box_impur, with rosenbluth trials.
  !! A series of positions for the non-impurity is generated, each
  !! deviating from the original position of the impurity wihin a radius
  !! of \ref input_parameters::rp_impur_rosbth. Each position has a weight
  !! \(w= \exp(-\beta E(trial)-E(old)) \). The finally attempted position is
  !! selected at random according to its weight. The introduced bias is corrected
  !! by computing the total weight of the old position and the new series of positions,
  !! these quantities modify the acceptance factor accordingly.
   implicit none
   type(PcmcState) :: st
   type(NrgRoutines) :: nr
   type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2,ener_dump
   logical :: bool,authtoselfam
   integer :: ip,ip1,ibox1,ibox2,i,idump,idump1,k,k_selected,counter=0
   integer :: fam1,fam2, ntrials_rosbth
   real(8) :: dump
   real(8),allocatable :: trial_pos(:,:),trial_prob(:)
   real(8) :: deltabmu, rosbth_fact,weight_new,weight_old,dumppos(3)
   type(NrgRec), allocatable :: Trial_ener(:)
   type(Particule) :: part1_0,part2_0
   save ::  trial_pos,trial_prob, Trial_ener
   counter=counter + 1
   ntrials_rosbth = st%inp%ntrials_rosbth
   if(.not.allocated(trial_pos)) then
      allocate(trial_pos(3,st%inp%ntrials_rosbth),trial_prob(ntrials_rosbth),Trial_ener(ntrials_rosbth))
   end if
   
   ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
   associate(labox2 => st%boxes(ibox2))
   ip1 = 1 + floor(ran2(st%rng)*labox2%n)
   fam2 = labox2%parts(ip1)%famille
   if(all(labox2%fam(2:st%dist%nfam)==0)) then
      echang2p_box_impur_rosbth=mv_cir_rej
      return
   end if
   do while(fam2==1)
      ip1 = 1 + floor(ran2(st%rng)*labox2%n)
      fam2 = labox2%parts(ip1)%famille
   end do
   
   ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
   if(st%ntotbox>1) then
      do while(ibox2 == ibox1)
         ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
      end do
   else
      echang2p_box_impur_rosbth=mv_cir_rej
      return
   end if
   associate(labox1 => st%boxes(ibox1))
   if(labox1%n<1.or.labox1%n<1) then
      echang2p_box_impur_rosbth=mv_cir_rej
      return
   end if
   
   ip = 1 + floor(ran2(st%rng)*labox1%n)
   if(fam2/=1 .and. labox1%fam(1)>0) authtoselfam = .true.        ! impurs are supposed to be fam 1
   if(authtoselfam) then
      do while(labox1%Parts(ip)%Famille /= 1 )
         ip = 1 + floor(ran2(st%rng)*labox1%n)
      enddo
   else
      echang2p_box_impur_rosbth=mv_cir_rej
      return
   end if
   
   ener1_1 = labox1%ener
   ener1_2 = labox2%ener
   ener0_1=ener1_1
   ener0_2=ener1_2
  
   associate(part1 => labox1%parts(ip), part2 => labox2%parts(ip1))
   ! part1 is impur
   fam1 = labox1%parts(ip)%famille
    
   part1_0 = part1
   part2_0 = part2

   ener1_1 = ener1_1 - nr%energy_1p(st,labox1,part1)
   ener_dump = nr%energy_1p(st,labox2,part2)
   ener1_2 = ener1_2 - ener_dump

   ! Old trials
    weight_old = exp(-ener_dump%ep/kT)
    do k=2,ntrials_rosbth   ! specific variable ?
       do i=1,3
          dump = part2_0%pos(i) + (st%inp%rp_impur_rosbth*(2.D0*ran2(st%rng)-1.D0))/labox2%tens(i,i)
          dumppos(i) = modulo(dump,1.D0)
          if(dumppos(i) >1.D0 .or. dumppos(i)<0.D0) then
             write(*,*) dump, dumppos(i) ; stop 666
          end if
       end do
       part2%pos(:) = dumppos(:)
       if(st%celllist) call celldec_update(part2,st%clistdecomp(ibox2))
       ener_dump = nr%energy_1p(st,labox2,part2)
       if(ener_dump%noverlap) then
          dump = exp(-ener_dump%ep/kT)
       else
          dump = 0.D0
       endif
       weight_old = weight_old +  dump
    end do
    part2%pos = part2_0%pos
    if(st%celllist) call celldec_update(part2,st%clistdecomp(ibox2))
    
    !Swap now...
    part1%numero = part2_0%numero
    part1%rayon = part2_0%rayon
    part1%famille = part2_0%famille
    part2%numero = part1_0%numero          ! part2 is impur (trial)
    part2%rayon = part1_0%rayon
    part2%famille = part1_0%famille
    ! "New" trials
    weight_new = 0.D0
    do k=1,ntrials_rosbth
       do i=1,3
          dump = part1_0%pos(i) + (st%inp%rp_IMPUR_ROSBTH*(2.D0*ran2(st%rng)-1.D0))/labox1%tens(i,i)
          trial_pos(i,k) = modulo(dump,1.D0)
       end do
       part1%pos(:) = trial_pos(1:3,k)
       if(st%celllist) call celldec_update(part1,st%clistdecomp(ibox1))
       ener_dump = nr%energy_1p(st,labox1,part1)
       if(ener_dump%noverlap) then
          trial_prob(k) = exp(-ener_dump%ep/kT)
          trial_ener(k) = ener_dump
       else
          trial_prob(k) = 0.D0
       endif
       weight_new = weight_new +  trial_prob(k)
    end do
 
    if(weight_new < 1.D-13) then ! Arbitrary thing Warning
       call restorepart()
       if(weight_new==0.D0) then
          echang2p_box_impur_rosbth=mv_hs_rej
       else
          echang2p_box_impur_rosbth=mv_nrj_rej
       end if
       return
    end if
    do k=2,ntrials_rosbth
       trial_prob(k)= trial_prob(k)+trial_prob(k-1)
    end do
    do k=1,ntrials_rosbth
       trial_prob(k) = trial_prob(k)/weight_new
    end do
    dump = ran2(st%rng)
    k_selected = ntrials_rosbth
    do k=1,ntrials_rosbth
       if(dump<=trial_prob(k)) then
          k_selected = k
          exit
       end if
    end do
    if(k_selected < 1 .or. k_selected > ntrials_rosbth) then
       write(*,*) real(trial_prob)
       stop 666
    end if
    part1%pos = trial_pos(1:3,k_selected)
    if(st%celllist) call celldec_update(part1,st%clistdecomp(ibox1))
    ener1_1%ep = ener1_1%ep + Trial_ener(k_selected)%ep
    rosbth_fact = weight_new/weight_old
    
    
    idump=max(labox1%fam(fam1),1)
    idump1=max(labox2%fam(fam2),1)
    deltabmu = log(((1.D0*(labox1%fam(fam2)+1))*(1.D0*(labox2%fam(fam1)+1)))/&
       ((1.D0*idump)*(1.D0*idump1)))
  
    if(trial_ener(k_selected)%noverlap) then
       ener1_2 = ener1_2 + nr%energy_1p(st,labox2,part1)
       if(ener1_2%noverlap) then
          if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)/kT + deltabmu - log(rosbth_fact),st%rng)) then
             labox1%ener=ener1_1
             labox2%ener=ener1_2
             labox1%fam(fam1)=labox1%fam(fam1)-1
             labox1%fam(fam2)=labox1%fam(fam2)+1
             labox2%fam(fam2)=labox2%fam(fam2)-1
             labox2%fam(fam1)=labox2%fam(fam1)+1
             echang2p_box_impur_rosbth = mv_acc
          else
             call restorepart()
             echang2p_box_impur_rosbth = mv_nrj_rej
          end if
       else
          call restorepart()
          echang2p_box_impur_rosbth = mv_hs_rej
       end if
    else
       call restorepart()
       echang2p_box_impur_rosbth = mv_hs_rej
    end if

    end associate
    end associate
    end associate

  contains

    subroutine restorepart()
      associate(part1 => st%boxes(ibox1)%parts(ip), part2 => st%boxes(ibox2)%parts(ip1))

      part1%numero = part1_0%numero
      part1%rayon = part1_0%rayon
      part1%famille = part1_0%famille
      part1%pos = part1_0%pos
      part2%numero = part2_0%numero          ! part2 is impur (trial)
      part2%rayon = part2_0%rayon
      part2%famille = part2_0%famille
      part2%pos = part2_0%pos

      if(st%celllist) call celldec_update(part2,st%clistdecomp(ibox2))
      if(st%celllist) call celldec_update(part1,st%clistdecomp(ibox1))

      end associate
    end subroutine restorepart
    
end function echang2p_box_impur_rosbth

integer function echang2p_box_impur_rosbth2(st, nr)
  !! Attempts to swap a particle with an "impurity"
  !! (ie the first particle family of particle, supposedly very small)
  !! between two different boxes, with rosenbluth trials (second method).
  !!
  !! Same principle that in \ref echang2p_box_impur_rosbth
  !! but in addition neighboring particles are virtually taken out the box
  !! and replaced into it at random around the new position of the swapped
  !! particle.
   use misc, only: search_neighbors
   implicit none
   type(PcmcState) :: st
   type(NrgRoutines) :: nr
   type(NrgRec) :: ener1_1,ener1_2,ener0_1,ener0_2,ener_dump
   logical :: bool,authtoselfam
   integer :: ip,ip1,ibox1,ibox2,i,idump,idump1,k,k_selected,counter=0
   integer :: fam1,fam2,j, ntrials_rosbth
   real(8) :: dump
   real(8),allocatable :: trial_pos(:,:),trial_prob(:)
   real(8) :: deltabmu, rosbth_fact,weight_new,weight_old,dumppos(3),weight_old_int,weight_new_int
   type(NrgRec), allocatable :: trial_ener(:)
   type(Particule) :: part1_0,part2_0
   integer,parameter :: nvoismax=30
   integer :: liste_voisins(nvoismax),nvois,nvoistrait   
   type(Particule) :: part_vois(nvoismax)
   save ::  trial_pos,trial_prob, trial_ener
   nvoistrait=0
   counter=counter + 1
   ntrials_rosbth = st%inp%ntrials_rosbth
   if(.not.allocated(trial_pos)) then
      allocate(trial_pos(3,st%inp%ntrials_rosbth),trial_prob(ntrials_rosbth),Trial_ener(st%inp%ntrials_rosbth))
   end if
   
   ibox2 = 1 + floor(ran2(st%rng)*st%ntotbox)
   associate(labox2 => st%boxes(ibox2))
   ip1 = 1 + floor(ran2(st%rng)*labox2%n)
   fam2 = labox2%parts(ip1)%famille
   if(all(labox2%fam(2:st%dist%nfam)==0)) then
      echang2p_box_impur_rosbth2=mv_cir_rej
      return
   end if
   do while(fam2==1)
      ip1 = 1 + floor(ran2(st%rng)*labox2%n)
      fam2 = labox2%parts(ip1)%famille
   end do
   
   ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
   if(st%ntotbox>1) then
      do while(ibox2 == ibox1)
         ibox1 = 1 + floor(ran2(st%rng)*st%ntotbox)
      end do
   else
      echang2p_box_impur_rosbth2=mv_cir_rej
      return
   end if
   associate(labox1 => st%boxes(ibox1))
   if(labox1%n<1.or.labox1%n<1) then
      echang2p_box_impur_rosbth2=mv_cir_rej
      return
   end if
   
   ip = 1 + floor(ran2(st%rng)*labox1%n)
   if(fam2/=1 .and. labox1%fam(1)>0) authtoselfam = .true.        !
   !impurs are supposed to be fam 1
   if(authtoselfam) then
      do while(labox1%Parts(ip)%Famille /= 1 )
         ip = 1 + floor(ran2(st%rng)*labox1%n)
      enddo
   else
      echang2p_box_impur_rosbth2=mv_cir_rej
      return
   end if
   
   ener1_1 = labox1%ener
   ener1_2 = labox2%ener
   ener0_1=ener1_1
   ener0_2=ener1_2
  
   associate(part1 => labox1%parts(ip), part2 => labox2%parts(ip1))
   !! part1 is impur
   fam1 = labox1%parts(ip)%famille
   
   part1_0 = part1
   part2_0 = part2

   ener_dump%ep=0.D0
   ener1_1 = ener1_1 - nr%energy_1p(st,labox1,part1)
   ener_dump = nr%energy_1p(st,labox1,part2)
   ener1_2 = ener1_2 - ener_dump

   ! Must do it before swap
   call search_neighbors(st,st%inp%rvois_impur_rosbth,part2,ibox2&
        &,nvois,liste_voisins,nvoismax)
   
   ! Old trials
    weight_old = exp(-ener_dump%ep/kT)
    do k=2,ntrials_rosbth   ! specific variable ?
       do i=1,3
          dump = part2_0%pos(i) + (st%inp%rp_impur_rosbth*(2.D0&
               &*ran2(st%rng)-1.D0))/labox2%tens(i,i)
          dumppos(i) = modulo(dump,1.D0)
          if(dumppos(i) >1.D0 .or. dumppos(i)<0.D0) then
              write(*,*) dump, dumppos(i) ; stop 666
          end if
       end do
       part2%pos(:) = dumppos(:)
       if(st%celllist) call celldec_update(part2,st&
            &%clistdecomp(ibox2))
       ener_dump = nr%energy_1p(st,labox2,part2)
       if(ener_dump%noverlap) then
          dump = exp(-ener_dump%ep/kT)
       else
          dump = 0.D0
       endif
       weight_old = weight_old +  dump
    end do
    part2%pos = part2_0%pos
    if(st%celllist) call celldec_update(part2,st%clistdecomp(ibox2))
    
    !Swap now...
    part1%numero = part2_0%numero
    part1%rayon = part2_0%rayon
    part1%famille = part2_0%famille
    part2%numero = part1_0%numero          ! part2 is impur (trial)
    part2%rayon = part1_0%rayon
    part2%famille = part1_0%famille
    ! "New" trials
    weight_new = 0.D0
    do k=1,ntrials_rosbth
       do i=1,3
          dump = part1_0%pos(i) + (st%inp%rp_IMPUR_ROSBTH*(2.D0&
               &*ran2(st%rng)-1.D0))/labox1%tens(i,i)
          trial_pos(i,k) = modulo(dump,1.D0)
       end do
       part1%pos(:) = trial_pos(1:3,k)
       if(st%celllist) call celldec_update(part1,st&
            &%clistdecomp(ibox1))
       ener_dump = nr%energy_1p(st,labox1,part1)
       if(ener_dump%noverlap) then
          trial_prob(k) = exp(-ener_dump%ep/kT)
          trial_ener(k) = ener_dump
       else
          trial_prob(k) = 0.D0
       endif
       weight_new = weight_new +  trial_prob(k)
    end do
 
    if(weight_new < 1.D-13) then ! Arbitrary thing Warning
       call RESTOREPART()
       if(weight_new==0.D0) then
          echang2p_box_impur_rosbth2=mv_hs_rej
       else
          echang2p_box_impur_rosbth2=mv_nrj_rej
       end if
       return
    end if
    do k=2,ntrials_rosbth
       trial_prob(k)= trial_prob(k)+trial_prob(k-1)
    end do
    do k=1,ntrials_rosbth
       trial_prob(k) = trial_prob(k)/weight_new
    end do
    dump = ran2(st%rng)
    k_selected = ntrials_rosbth
    do k=1,ntrials_rosbth
       if(dump<=trial_prob(k)) then
          k_selected = k
          exit
       end if
    end do
    if(k_selected < 1 .or. k_selected > ntrials_rosbth) then
       write(*,*) real(trial_prob)
       stop 666
    end if
    part1%pos = trial_pos(1:3,k_selected)
    if(st%celllist) call celldec_update(part1,st%clistdecomp(ibox1))
    ener1_1%ep = ener1_1%ep + Trial_ener(k_selected)%ep

    ! Biased displacements of surrounding parts of part1 ( which is
    ! now part2 )
    ! Old trials
    do i=1,nvois
       part_vois(i) = labox2%parts(liste_voisins(i))
       associate(part_temp => labox2%parts(liste_voisins(i)))

       ener_dump = nr%energy_1p(st,labox2,part_temp)
       
       weight_old_int = exp(-ener_dump%ep/kT)
       do k=2,ntrials_rosbth   ! specific variable ?
          do j=1,3
             dump = part_vois(i)%pos(i) + (st%inp%rp_impur_rosbth&
                  &*(2.D0*ran2(st%rng)-1.D0))/labox2%tens(i,i)
             dumppos(i) = modulo(dump,1.D0)
          end do
          part_temp%pos(:) = dumppos(:)
          if(st%celllist) call celldec_update(part_temp,st&
               &%clistdecomp(ibox2))
          ener_dump = nr%energy_1p(st,labox2,part_temp)
          if(ener_dump%noverlap) then
             dump = exp(-ener_dump%ep/kT)
          else
             dump = 0.d0
          endif
          weight_old_int = weight_old_int +  dump
       end do
       weight_old = weight_old*weight_old_int
       end associate
    end do

    call search_neighbors(st,st%inp%rvois_impur_rosbth,part1,ibox1&
         &,nvois,liste_voisins,nvoismax)
    ! New trials
    nvoistrait=0
    do i=1,nvois
       nvoistrait=nvoistrait+1
       part_vois(i) = labox1%parts(liste_voisins(i))
       associate(part_temp => labox1%parts(liste_voisins(i)))

       ener1_1 = ener1_1 - nr%energy_1p(st,labox1,part_temp)
       
       weight_new_int = 0.D0
       do k=1,ntrials_rosbth
          do j=1,3
             dump = part_vois(j)%pos(j) + (st%inp%rp_impur_rosbth&
                  &*(2.d0*ran2(st%rng)-1.d0))/labox1%tens(i,i)
             trial_pos(j,k) = modulo(dump,1.d0)
          end do
          part_temp%pos(:) = trial_pos(1:3,k)
          if(st%celllist) call celldec_update(part_temp,st&
               &%clistdecomp(ibox1))
          ener_dump = nr%energy_1p(st,labox1,part_temp)
          if(ener_dump%noverlap) then
             trial_prob(k) = exp(-ener_dump%ep/kT)
             trial_ener(k) = ener_dump
          else
             trial_prob(k) = 0.D0
          endif
          weight_new_int = weight_new_int +  trial_prob(k)
       end do
 
       if(weight_new_int < 1.d-13) then ! Arbitrary thing Warning
          call restorepart()
          if(weight_new_int==0.D0) then
             echang2p_box_impur_rosbth2=mv_hs_rej
          else
             echang2p_box_impur_rosbth2=mv_nrj_rej
          end if
          return
       end if
       do k=2,ntrials_rosbth
          trial_prob(k)= trial_prob(k)+trial_prob(k-1)
       end do
       do k=1,ntrials_rosbth
          trial_prob(k) = trial_prob(k)/weight_new
       end do
       dump = ran2(st%rng)
       k_selected = ntrials_rosbth
       do k=1,ntrials_rosbth
          if(dump<=trial_prob(k)) then
             k_selected = k
             exit
          end if
       end do
       if(k_selected < 1 .or. k_selected > ntrials_rosbth) then
          write(*,*) real(trial_prob)
          stop 666
       end if
       part_temp%pos = trial_pos(1:3,k_selected)
       if(st%celllist) call celldec_update(part_temp,st&
            &%clistdecomp(ibox1))
       ener1_1 = ener1_1 + trial_ener(k_selected)

       weight_new = weight_new * weight_new_int
       end associate
    end do
    
    rosbth_fact = weight_new/weight_old
    idump=max(labox1%fam(fam1),1)
    idump1=max(labox2%fam(fam2),1)
    deltabmu = log(((1.d0*(labox1%fam(fam2)+1))*(1.d0*(labox2&
         &%fam(fam1)+1)))/ ((1.d0*idump)*(1.d0*idump1)))
  
    if(ener1_1%noverlap) then
       ener1_2 = ener1_2 + nr%energy_1p(st,labox2,part1)
       if(ener1_2%noverlap) then
          if(metropolis((ener1_1%ep+ener1_2%ep-ener0_1%ep-ener0_2%ep)&
               &/kT + deltabmu- log(rosbth_fact),st%rng)) then
             labox1%ener=ener1_1
             labox2%ener=ener1_2
             labox1%fam(fam1)=labox1%fam(fam1)-1
             labox1%fam(fam2)=labox1%fam(fam2)+1
             labox2%fam(fam2)=labox2%fam(fam2)-1
             labox2%fam(fam1)=labox2%fam(fam1)+1
             echang2p_box_impur_rosbth2 = mv_acc
          else
             call restorepart()
             echang2p_box_impur_rosbth2 = mv_nrj_rej
          end if
       else
          call restorepart()
          echang2p_box_impur_rosbth2 = mv_hs_rej
       end if
    else
       call restorepart()
       echang2p_box_impur_rosbth2 = mv_hs_rej
    end if

    end associate
    end associate
    end associate

  contains

    subroutine restorepart()
      associate(labox1 => st%boxes(ibox1), labox2 => st%boxes(ibox2), &
        part1 => st%boxes(ibox1)%parts(ip), part2 => st%boxes(ibox2)%parts(ip1))

      part1%numero = part1_0%numero
      part1%rayon = part1_0%rayon
      part1%famille = part1_0%famille
      part1%pos = part1_0%pos
      part2%numero = part2_0%numero          ! part2 is impur (trial)
      part2%rayon = part2_0%rayon
      part2%famille = part2_0%famille
      part2%pos = part2_0%pos
      if(st%celllist) call celldec_update(part2,st&
           &%clistdecomp(ibox2))
      if(st%celllist) call celldec_update(part1,st&
           &%clistdecomp(ibox1))

      do i=1,nvoistrait
         labox1%parts(liste_voisins(i))%numero = part_vois(i)%numero
         labox1%parts(liste_voisins(i))%rayon = part_vois(i)%rayon
         labox1%parts(liste_voisins(i))%famille = part_vois(i)%famille
         labox1%parts(liste_voisins(i))%pos = part_vois(i)%pos
         if(st%celllist) call celldec_update(labox1&
              &%parts(liste_voisins(i)),st%clistdecomp(ibox1))
      end do
      
      end associate
    end subroutine restorepart
    
end function echang2p_box_impur_rosbth2 

end module spemoves

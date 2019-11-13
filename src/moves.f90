module moves
  !! Contains helpers for user-defined moves
  !!
  !! Before running the simulation, call [[add_move]] to register your moves.
  !! The sum of all probabilities must be inferior to 1. The remaining probability is
  !! the probability of doing a simple particle translation.
  !!
  !! The move procedure should at least return [[moves:mv_acc]] when the move is accepted
  !! and something else otherwise.
  !!
  !!### Example
  !!
  !!```fortran
  !!use moves
  !!integer function cluster_move(st)
  !!  type(PcmcState), intent(in) :: st
  !!  ! implementation left to the reader
  !!  cluster_move = MV_ACC
  !!end function
  !!
  !!type(mcsimulation) :: simu
  !!call add_move(simu%moves,"cluster move", 0.1d0, cluster_move)
  !!```

  use state, only: PcmcState
  use nrg
  implicit none
  
  interface
    integer function handler_type(st,nr)
      import
      type(PcmcState) :: st
      type(NrgRoutines) :: nr
    end function handler_type
  end interface

  type Mover
    !! Contains info, statistic and handler on MC moves
    integer(8) :: comp !! number of  attempted moves
    integer(8) :: acc    !! number of  accepted moves
    integer(8) :: hcrej  !! number of  rejected moves (hard overlap)
    integer(8) :: cirrej  !! number of rejected moves (neither for overlap nor for energy)
    character(48) :: name !! name of the move
    procedure(handler_type),pointer,nopass :: handler  !! function that attempts the move
  end type Mover

  integer, parameter, private ::  min_mouvs=15

  type MCMoves
    !! Set of moves
    integer :: nmoves = 0
    real(8), allocatable :: prob_mouvs(:)
    type(Mover), allocatable :: mv(:)
  end type MCMoves

  ! Conventions
  integer,parameter :: mv_acc=0      !! Accepted move
  integer,parameter :: mv_hs_rej=1   !! Move rejected for hard sphere overlap
  integer,parameter :: mv_nrj_rej=2  !! Move rejected for energy condition
  integer,parameter :: mv_bc_rej=4   !! Move rejected for out of bondaries condition
  integer,parameter :: mv_cir_rej=3  !! Move rejected for misc circumstances

contains

  subroutine deallocate_moves(mvs)
    type(MCMoves) :: mvs
    if(allocated(mvs%mv)) then
      deallocate(mvs%mv, mvs%prob_mouvs)
    end if
  end subroutine

  subroutine moves_reset(mvs)
    !! Reset the number of additional moves to 0.
    type(MCMoves) :: mvs
    mvs%nmoves = 0
  end subroutine moves_reset

  subroutine add_move(mvs,name,prob,move_proc)
    !! Registers a move in `mvs`.
    type(MCMoves) :: mvs !! container of moves
    real(8), intent(in) :: prob !! Probability of the move (in [0,1[)
    character(*) :: name      !! Name of the move (used for printing stats)
    procedure(handler_type) :: move_proc !! procedure that performs the move itself (signature integer function(PcmcState))
    integer :: nmoves, mvcap, omvcap
    type(Mover), allocatable :: temp_mv(:)
    real(8), allocatable :: temp_pm(:)
    if(.not. allocated(mvs%mv)) then
      allocate(mvs%mv(min_mouvs), mvs%prob_mouvs(min_mouvs))
    end if
    if(prob > 0) then
      nmoves = mvs%nmoves+1

      omvcap = ubound(mvs%mv,1)
      if(nmoves > omvcap) then
        call move_alloc(mvs%mv, temp_mv)
        call move_alloc(mvs%prob_mouvs, temp_pm)
        mvcap = 1 + floor(omvcap*1.5)
        allocate(mvs%mv(mvcap), mvs%prob_mouvs(mvcap))
        mvs%mv(1:omvcap) = temp_mv(1:omvcap)
        mvs%prob_mouvs(1:omvcap) = temp_pm(1:omvcap)
        deallocate(temp_mv, temp_pm)
      end if

      if(nmoves == 1) then
        mvs%prob_mouvs(1) = 0.0 ! 1st move is default move 
      else
        mvs%prob_mouvs(nmoves) = mvs%prob_mouvs(nmoves-1) + prob
      end if

      mvs%mv(nmoves)%name = name
      mvs%mv(nmoves)%comp=0
      mvs%mv(nmoves)%acc=0
      mvs%mv(nmoves)%hcrej=0
      mvs%mv(nmoves)%cirrej=0

      mvs%mv(nmoves)%handler => move_proc
      mvs%nmoves = nmoves
    end if
  end subroutine add_move

  subroutine reset_all_moves(mvs)
    !! Reset `mvs`. Useful to make a new run.
    type(MCMoves) :: mvs
    integer :: i
    do i=1,mvs%nmoves
      mvs%mv(i)%comp=0
      mvs%mv(i)%acc=0
      mvs%mv(i)%hcrej=0
      mvs%mv(i)%cirrej=0
    end do
  end subroutine reset_all_moves

  subroutine new_movinst(st,nr,mouv)
    !! Make a move according to `mouv`
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    type(Mover) :: mouv
    integer :: stat

    stat = mouv%handler(st,nr)
    mouv%comp = mouv%comp + 1
    if(stat == mv_acc) then
      mouv%acc = mouv%acc + 1
    elseif(stat == mv_hs_rej) then
      mouv%hcrej = mouv%hcrej + 1
    elseif(stat == mv_cir_rej) then
      mouv%cirrej = mouv%cirrej + 1       
    end if
  end subroutine new_movinst

  integer function select_move(mvs,pppp)
    !! Selects a move with random number `pppp`. (helper function).
    !! Returns index of the selected move.
    type(MCMoves) :: mvs
    real(8),intent(in) :: pppp
    integer :: ii

    select_move=1 ! Defaults on first move
    do ii=2,mvs%nmoves
      if(pppp<mvs%prob_mouvs(ii)) then
        select_move=ii
        exit
      end if
    end do

  end function select_move

  subroutine change_position_random(x, relRay, outfo, st)
    !! Changes a position `x` at random (within an interval).
    !! `outfo` returns information about the move: 0 means a wall has been hit,
    !! 1 is for a regular translation, 2 is for a translation near a wall.
    use random2
    implicit none
    type(PcmcState) :: st
    real(8),intent(in) :: relRay
    integer, intent(out) :: outfo
    real(8),intent(inout) :: x(3)
    integer :: j
    real(8) :: dx, dump,dump2, deltawall
    if(st%mode_box == 2) then
      deltawall = st%rpar%dx_trans + relRay
      if(x(3) < deltawall .or. 1.d0-x(3) < deltawall) then
        ! Special treatment near a wall. relRay should be the (reduced) radius of the particle.
        dx = st%rpar%dx_trans_wall
        outfo = 2
      else
        dx = st%rpar%dx_trans
        outfo = 1
      end if
      !PBC on X,Y
      do j=1,2
        dump2=x(j)+dx*(ran2(st%rng)-0.5d0)
        dump=modulo(dump2,1.d0)
        x(j)=dump
      enddo
      ! Hard Z Wall
      dump2=x(3)+dx*(ran2(st%rng)-0.5d0)
      if(dump2 < 0.d0  .or. dump2 >  1.d0) then
        outfo = 0
        return
      endif
    else
      dx = st%rpar%dx_trans
      outfo = 1
      do j=1,3
        dump2=x(j)+dx*(ran2(st%rng)-0.5d0)
        dump=modulo(dump2,1.d0)
        x(j)=dump
      enddo
    end if
  
  end subroutine change_position_random

  integer function translation1P(st, nr)
    !! Attempts a simple translation of a particle in a box.
    use state
    implicit none
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    type(NrgRec) :: ener1,ener0
    integer :: ip,ibox,counter=0, outfo_rcp
    real(8) :: store_pos(3),pp(3)
    counter=counter+1
    ibox = 1 + floor(ran2(st%rng)*st%ntotbox)
    ip = 1 + floor(ran2(st%rng)*st%boxes(ibox)%n)
    associate(part => st%boxes(ibox)%Parts(ip))

    store_pos=part%pos
    
    ener0 = nr%energy_1p(st,st%boxes(ibox),part)

    pp = part%pos
    call change_position_random(pp, part%rayon/st%boxes(ibox)%tens(3,3), outfo_rcp,st)
    if(outfo_rcp==0) then
      translation1P=4
      return
    end if
    call change_part_pos(st,ibox, ip, pp, store_pos)
    ener1 = nr%energy_1p(st,st%boxes(ibox),part)
    if(ener1%noverlap) then
      if(metropolis((ener1%ep-ener0%ep)/kT,st%rng)) then
        st%boxes(ibox)%ener%ep = st%boxes(ibox)%ener%ep + ener1%ep - ener0%ep
        translation1P=0
      else
        call change_part_pos(st,ibox, ip, store_pos, pp)
        translation1P=2
      endif
    else
      call change_part_pos(st,ibox, ip, store_pos, pp)
      translation1P=1
    endif

    st%sdata%trans_eff = st%sdata%trans_eff + 1
    if(translation1P==0) st%sdata%reuss_trans=st%sdata%reuss_trans + 1
    if(st%mode_box == 2) then
      st%transss = st%transss + 1
      if(outfo_rcp == 2) then
        st%trans_wall = st%trans_wall + 1
        if(translation1P==0) st%reuss_trans_wall = st%reuss_trans_wall + 1
      else if(outfo_rcp == 1) then
        st%trans_outwall = st%trans_outwall + 1
        if(translation1P==0) st%reuss_trans_outwall = st%reuss_trans_outwall + 1
      end if
    end if

    end associate
  end function translation1P

  subroutine change_part_pos(st, box, ip, new_pos, old_pos)
    !! Changes the position of particle. Returns old position in `old_pos`.
    use celldec
    type(PcmcState) :: st
    integer,intent(in) :: ip,box
    real(8), intent(in) :: new_pos(3)
    real(8), intent(out) :: old_pos(3)
    old_pos = st%boxes(box)%parts(ip)%pos
    st%boxes(box)%parts(ip)%pos = new_pos
    if(st%celllist) call celldec_update(st%boxes(box)%parts(ip),st%clistdecomp(box))
  end subroutine change_part_pos

  subroutine change_part_radius(st, box, ip, new_radius)
    !! Changes the radius of particle `ip`.
    type(PcmcState) :: st
    integer,intent(in) :: ip, box
    real(8), intent(in) :: new_radius
    integer :: old_fam, fam, j
    real(8) :: mr
    st%boxes(box)%parts(ip)%rayon = new_radius
    old_fam = st%boxes(box)%parts(ip)%famille
    ! Updates family within the predefined distribution.
    ! ⚠ Supposes a sorted distribution, may not work.
    if(st%dist%nfam <= 1) then
      fam = 1
    else if(new_radius <= st%dist%rad(1)) then
      fam = 1
    else if(new_radius >= st%dist%rad(st%dist%nfam)) then
      fam = st%dist%nfam
    else
      do j=1, st%dist%nfam-1, 1
        if(st%dist%rad(j) <= new_radius .and. new_radius < st%dist%rad(j+1)) exit
      end do
      mr = 0.5d0*(st%dist%rad(j+1) + st%dist%rad(j))
      if(new_radius <= mr) then
        fam = j
      else
        fam = j+1
      end if
    end if
    if(old_fam /= fam) call change_part_family(st, box, ip, fam)
  end subroutine change_part_radius

  subroutine change_part_echarge(st, box, ip, new_charge)
    !! Changes the internal charge of particle `ip`.
    !! @Warning 
    !! For yukawa model it is (effective charge) × exp(κ R)/(1+κ R).
    !! In general adapt to the model used.
    !! @endwarning
    type(PcmcState) :: st
    integer,intent(in) :: ip, box
    real(8), intent(in) :: new_charge
    st%boxes(box)%parts(ip)%ech = new_charge
  end subroutine change_part_echarge
  
  subroutine change_part_family(st, box, ip, new_family)
    !! Changes the family of a particle `ip`.
    type(PcmcState) :: st
    integer,intent(in) :: ip, box
    integer, intent(in) :: new_family
    integer :: old_fam
    old_fam = st%boxes(box)%parts(ip)%famille
    ! should check if famille is in range
    if(old_fam /= new_family) then
      st%boxes(box)%parts(ip)%famille = new_family
      st%boxes(box)%fam(old_fam) = st%boxes(box)%fam(old_fam) - 1
      st%boxes(box)%fam(new_family) = st%boxes(box)%fam(new_family) + 1
    end if
  end subroutine change_part_family

  subroutine remove_part_box(st, box, ip)
    !! Removes particle `ip` from box `box`.
    use celldec
    implicit none
    type(PcmcState) :: st
    integer,intent(in) :: box, ip
    type(Particule) :: part
    call remove_part_box_nocelllist(st, box, ip, part)
    if(st%celllist) call celldec_update_all(st%clistdecomp(box), st%boxes(box)%n, st%boxes(box)%parts)
  end subroutine remove_part_box

  subroutine remove_part_box_nocelllist(st,box,ip, part)
    !> Removes particle *ip* from box *box*.
    !! @Warning 
    !! This version does not update the cell decomposition.
    !! @endwarning
    use celldec, only: Particule
    implicit none
    type(PcmcState) :: st
    integer, intent(in) :: ip,box
    type(Particule), intent(out) :: part
    integer :: i,f

    part = st%boxes(box)%parts(ip)
    f = st%boxes(box)%parts(ip)%famille

    st%boxes(box)%fam(f) = st%boxes(box)%fam(f) - 1
    st%boxes(box)%n = st%boxes(box)%n - 1
    ! Shuffle. Invalidates all pointers (and reference indices) to parts(ip:n) !
    ! And O(n) operation. Alternative is to make holes.
    do i=ip,st%boxes(box)%n
      st%boxes(box)%parts(i) = st%boxes(box)%parts(i+1)
    end do
  end subroutine remove_part_box_nocelllist

  subroutine add_part_box(st, box, part)
    !> Adds particle *part* to box *box*.
    use celldec
    implicit none
    type(PcmcState) :: st
    type(Particule),intent(in) :: part  ! particle
    integer,intent(in) :: box  !! box number
    call add_part_box_nocelllist(st, box, part)
    if(st%celllist) then
      call celldec_add(st%boxes(box)%parts(st%boxes(box)%n),st%clistdecomp(box))
    end if
  end subroutine add_part_box

  subroutine add_part_box_nocelllist(st,box,part)
    !> Adds particle *part* to box *box*.
    !! @Warning
    !! This version does not update the cell decomposition.
    !! @endwarning
    use celldec
    use iso_fortran_env
    implicit none
    type(PcmcState) :: st
    type(Particule),intent(in) :: part
    integer,intent(in) :: box

    integer :: no
    no = st%boxes(box)%n
    if((no+1) > ubound(st%boxes(box)%parts,1)) then
      write(error_unit, '(A,I0,A)') "Error: Reached the maximum number of particles (", ubound(st%boxes(box)%parts,1),")."
      write(error_unit,'(A)') "       You can increase 'max_number_of_particles' in your configuration."
      stop 9
    end if
    st%boxes(box)%n=no+1
    st%boxes(box)%parts(no+1) = part
    st%boxes(box)%fam(part%famille) = st%boxes(box)%fam(part%famille) + 1
  end subroutine add_part_box_nocelllist

  subroutine swap_2parts(st,ibox, part1, part2)
    !! Swaps two particles.
    use celldec
    type(PcmcState) :: st
    integer, intent(in) :: ibox  !! box number
    type(Particule), intent(inout) :: part1, part2
    real(8) :: store_pos(3)
    store_pos=part1%pos
    part1%pos=part2%pos
    part2%pos=store_pos
    if(st%celllist) call celldec_echange(part1,part2,st%clistdecomp(ibox))
  end subroutine swap_2parts

  subroutine swap_2parts_2box(st,ibox1, ibox2, part1, part2)
    !! Swaps two particles between two boxes.
    use celldec
    type(PcmcState) :: st
    integer, intent(in) :: ibox1, ibox2
    type(Particule), intent(inout) :: part1, part2
    type(Particule) :: store_part
    integer :: fam1, fam2
    fam1 = part1%famille
    fam2 = part2%famille
    store_part = part1
    ! Note: Reference to particles (index or pointer) would points to a different particle
    call part_transfer_props(part2, part1)
    call part_transfer_props(store_part, part2)
    ! Needs to update box distributions
    st%boxes(ibox1)%fam(fam1)=st%boxes(ibox1)%fam(fam1)-1
    st%boxes(ibox1)%fam(fam2)=st%boxes(ibox1)%fam(fam2)+1
    st%boxes(ibox2)%fam(fam2)=st%boxes(ibox2)%fam(fam2)-1
    st%boxes(ibox2)%fam(fam1)=st%boxes(ibox2)%fam(fam1)+1
  end subroutine swap_2parts_2box

  subroutine transfer_part_b2b(st,part,box1,box2)
    !! Transfer a particle from box `box1` to box `box2`.
    use celldec
    use iso_fortran_env
    implicit none
    type(PcmcState),target :: st
    type(Particule),target :: part
    integer, intent(in) :: box1,box2

    type(Particule),pointer :: partN
    integer :: fam

    ! Add part to box2
    fam = part%famille
    st%boxes(box2)%n = st%boxes(box2)%n + 1
    if(st%boxes(box2)%n > ubound(st%boxes(box2)%parts,1)) then
      write(error_unit, '(A,I0,A)') "Error: Reached the maximum number of particles (", ubound(st%boxes(box2)%parts,1),")."
      write(error_unit,'(A)') "       To work around this you can increase 'max_number_of_particles' in your configuration."
      stop 3
    end if
    st%boxes(box2)%fam(fam) = st%boxes(box2)%fam(fam) + 1
    call part_transfer_props(part, st%boxes(box2)%parts(st%boxes(box2)%n))
    st%boxes(box2)%parts(st%boxes(box2)%n)%pos = part%pos
    if(st%celllist) call celldec_add(st%boxes(box2)%parts(st%boxes(box2)%n),st%clistdecomp(box2))

    ! Then the trick is to transfer the last particle of box1 to part location
    partN => st%boxes(box1)%parts(st%boxes(box1)%n)
    if(.not. associated(partN, part)) then
      call part_transfer_props(partN, part)
      part%pos = partN%pos
    end if
    ! Then the new last particle can be "removed".
    st%boxes(box1)%n = st%boxes(box1)%n - 1
    st%boxes(box1)%fam(fam) = st%boxes(box1)%fam(fam) - 1
    ! Updates cell decomposition
    if(st%celllist) then
      if(associated(part%precedent)) then
        part%precedent%suivant => part%suivant
      else
        st%clistdecomp(box1)%cells(part%cellule)%hoc => part%suivant
      end if
      if(associated(part%suivant)) part%suivant%precedent => part%precedent

      if(.not.associated(partN,part)) then
        part%precedent => partN%precedent
        part%suivant => partN%suivant
        if(associated(partN%suivant)) partN%suivant%precedent => part
        if(associated(partN%precedent)) then
          partN%precedent%suivant => part
        else
          st%clistdecomp(box1)%cells(partN%cellule)%hoc => part
        end if
        part%cellule = partN%cellule
      end if
    end if

  end subroutine transfer_part_b2b

end module moves

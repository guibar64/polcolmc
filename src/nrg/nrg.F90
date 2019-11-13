
module nrg
!! Sets up builtin energy routines
!!
use state, only: PcmcState, NrgRec, Box, Pressure
use celldec, only: Particule
use none, only: none_initpot => initpot, none_nrg_1p => nrg_1p, &
     none_nrg_box => nrg_box, none_potent => potential, none_force => force 
use yukawa, only: yukawa_initpot => initpot, yukawa_nrg_1p => nrg_1p, &
     yukawa_nrg_box => nrg_box, yukawa_potent => potential, yukawa_force => force
use yukawa2, only: yukawa2_initpot => initpot, yukawa2_nrg_1p => nrg_1p, &
     yukawa2_nrg_box => nrg_box, yukawa2_potent => potential, yukawa2_force => force
use fennel, only: fennel_initpot => initpot, fennel_nrg_1p => nrg_1p, &
     fennel_nrg_box => nrg_box, fennel_potent => potential, fennel_force => force

implicit none
abstract interface
  real(8) function potent_proc(r,i,j,st)
    import
    integer,intent(in) :: i,j
    real(8),intent(in) :: r
    type(PcmcState) :: st
  end function
  real(8) function force_proc(r,i,j,st)
    import
    integer,intent(in) :: i,j
    real(8),intent(in) :: r
    type(PcmcState) :: st
  end function
  function nrg_1p_proc(st,b,par)
    import
    type(NrgRec) :: nrg_1p_proc
    type(PcmcState), intent(in) :: st
    type(Box),intent(in),target :: b   !! simulation box
    type(Particule), intent(in) :: par !! the particle
  end function
  function nrg_box_proc(st,b,press)
    import
    type(NrgRec) :: nrg_box_proc
    type(PcmcState), intent(in) :: st
    type(Box), intent(in) :: b                     !! simulation box
    type(Pressure),intent(out), optional :: press  !! Pressure
  end function
  subroutine initpot_proc(st)
    import
    type(PcmcState) :: st
  end subroutine
end interface

type NrgRoutines
  procedure(potent_proc), nopass, pointer :: potential       !! pair potential 
  procedure(force_proc), nopass, pointer :: force            !! pair force
  procedure(nrg_1p_proc), nopass, pointer :: energy_1p       !! energy of one particle
  procedure(nrg_box_proc), nopass, pointer :: energy_box     !! energy of a simulation box
  procedure(initpot_proc), nopass, pointer :: init_potential !! initialize potential
end type

private :: eff_charge

contains

subroutine modes_affect_pot(nr, s_mod, m)
  !! Affects potential modes from ``s_mod``.
  use iso_fortran_env

  implicit none
  type(NrgRoutines) :: nr
  character(*),intent(in) :: s_mod
  integer, intent(out) :: m
  select case(s_mod)
  case("none")
    m = 0
    call set_nrg_routines(nr, none_potent, none_force, none_nrg_1p, none_nrg_box, none_initpot)
  case("yukawa")
    m = 1
    call set_nrg_routines(nr, yukawa_potent, yukawa_force, yukawa_nrg_1p, yukawa_nrg_box, yukawa_initpot)
  case("yukawa2")
    m = 2
    call set_nrg_routines(nr, yukawa2_potent, yukawa2_force, yukawa2_nrg_1p, yukawa2_nrg_box, yukawa2_initpot)
  case("fennel")
    m = 3
    call set_nrg_routines(nr, fennel_potent, fennel_force, fennel_nrg_1p, fennel_nrg_box, fennel_initpot)
  case("user")
    m = 4
    if(.not. associated(nr%potential)) then
       write(error_unit,'(A)') "Selected user-defined potential but routines do not seem associated. &
&Have you called 'set_nrg_routines' before initializing the simulation ?"
       stop 2
     end if
  case default
     write(error_unit,'(A,A,A)') "Error: Unkown potential '", s_mod, "'."
     stop 2
  end select
end subroutine modes_affect_pot

subroutine set_nrg_routines(nr, potent_routine, force_routine, nrg_1p_routine, nrg_box_routine, init_potential_routine)
  !! Set custom routines that computes energy / pressures. 
  !! @warning
  !! If `potential_type` is not set to "user", [[modes_affect_pot]] may override the routines
  !! @endwarning
  type(NrgRoutines) :: nr
  procedure(potent_proc) :: potent_routine
  procedure(force_proc) :: force_routine
  procedure(nrg_1p_proc) :: nrg_1p_routine
  procedure(nrg_box_proc) :: nrg_box_routine
  procedure(initpot_proc) :: init_potential_routine

  nr%potential => potent_routine
  nr%force =>  force_routine
  nr%energy_1p => nrg_1p_routine
  nr%energy_box => nrg_box_routine
  nr%init_potential => init_potential_routine

end subroutine set_nrg_routines


real(8) function eff_charge(z, kp, r, old)
  real(8), intent(in) :: z, kp, r
  logical :: old
  if(old) then
    eff_charge = z*&
            sinh(kp*r)/(kp*r)
  else
    eff_charge = z*&
            exp(kp*r)/(1+kp*r)
  end if
end function

subroutine nrg_init(st, nr)
  !! Initialize energy-related variables and calls
  !! potential-specific initialization.
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  integer :: mode,i
  real(8) :: kp

  call modes_affect_pot(nr, st%inp%smod_pot, st%mode_pot)

  mode = st%mode_box + st%mode_dec*100

  st%ndat%rcutoff_max = huge(1.d0)
  kp=1.d0/st%inp%debye_length
  allocate(st%ndat%eff_charge(st%dist%nfam))
  if(st%mode_pot == 1 .or. st%mode_pot == 2 ) then ! yukawa(2)
    do i=1,st%dist%nfam
      st%ndat%eff_charge(i) = eff_charge(st%dist%charge(i), kp, st%dist%rad(i), st%inp%old_fashion)
    end do
  end if
  st%ndat%inv_debye_length=kp
  st%ndat%elpot_pref=st%inp%length_prefactor
  st%mode_combo = mode

  call nr%init_potential(st)

end subroutine nrg_init

end module nrg

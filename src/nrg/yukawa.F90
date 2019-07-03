module yukawa
  !! Yukawa potential
  use state
  implicit none

  private :: potent_r, potent_r2, force_r

! The nrg_(1p|box)_* should be private to the module 
#include "nrgprivs.f90"

contains

subroutine initpot(st)
  type(PcmcState) :: st
end subroutine initpot


real(8) function potent_r(r, nd)
  implicit none
  real(8), intent(in) :: r
  type(NrgData), intent(in) :: nd

  potent_r=exp(-r*nd%inv_debye_length)/r
end function potent_r

real(8) function potent_r2(r2, nd)
  implicit none
  real(8), intent(in) :: r2
  type(NrgData), intent(in) :: nd
  real(8) :: r

  r=sqrt(r2)  
  potent_r2=exp(-r*nd%inv_debye_length)/r
end function potent_r2

real(8) function force_r(r, nd)
  implicit none
  real(8), intent(in) :: r
  type(NrgData), intent(in) :: nd

  real(8) :: dump
   

  dump=1/r 
  force_r=dump*(dump+nd%inv_debye_length)*exp(-r*nd%inv_debye_length)
end function force_r

real(8) function potential(r,i,j,st)
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: r
  type(PcmcState) :: st
  potential = st%ndat%elpot_pref * st%ndat%eff_charge(i)*st%ndat%eff_charge(j)*potent_r(r, st%ndat)
end function potential

real(8) function force(r,i,j,st)
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: r
  type(PcmcState) :: st
  force = st%ndat%elpot_pref * st%ndat%eff_charge(i)*st%ndat%eff_charge(j)*force_r(r, st%ndat)
end function force

! nrg_(1p|box)_* implementations that will use the potent_*, force_* previously defined 
#include "nrgimp.f90"

end module

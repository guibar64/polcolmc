!$$$$$$$$$$$$$$$$$ Fennel potential $$$$$$$$$$$$$$$$
!
!  u_ij = q_i * q_j * l_b * (1/r + r/r_c^2 - 2 / r_c)
!
module fennel
!! Fennel potential
use state
implicit none

private :: potent_r, potent_r2, force_r

! The nrg_(1p|box)_* should be private to the module
#include "nrgprivs.f90"

contains

subroutine initpot(st)
  use iso_fortran_env
  type(PcmcState) :: st
  integer :: i

  do i=1,st%dist%nfam
      st%ndat%eff_charge(i)=st%dist%charge(i)
  end do
  if(st%inp%rcutoff < 0) then
    write(error_unit, *) "Error: For potential 'fennel' you must specify a 'cutoff radius' > 0 in the configuration file."
    stop 6
  end if
  st%ndat%inv_rc2 = 1.d0/st%inp%rcutoff**2
  st%ndat%minus2_inv_rc = -2.d0/st%inp%rcutoff

end subroutine

real(8) function potent_r(r, nd)
  implicit none
  real(8), intent(in) :: r
  type(NrgData), intent(in) :: nd

  potent_r=(1._8/r + r*nd%inv_rc2 + nd%minus2_inv_rc)
end function

real(8) function potent_r2(r2, nd)
  implicit none
  real(8), intent(in) :: r2
  type(NrgData), intent(in) :: nd
  real(8) :: r

  r=sqrt(r2)
  potent_r2=(1/r + r*nd%inv_rc2 + nd%minus2_inv_rc)
end function


real(8) function force_r(r, nd)
  implicit none
  real(8), intent(in) :: r
  type(NrgData), intent(in) :: nd
  real(8) :: dump

  dump=1.d0/r**2
  force_r=(dump - nd%inv_rc2**2)
end function

real(8) function potential(r,i,j,st)
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: r
  type(PcmcState) :: st
  potential = st%ndat%elpot_pref*st%ndat%eff_charge(i)*st%ndat%eff_charge(j)*potent_r(r, st%ndat)
end function potential

real(8) function force(r,i,j,st)
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: r
  type(PcmcState) :: st
  force = st%ndat%elpot_pref*st%ndat%eff_charge(i)*st%ndat%eff_charge(j)*force_r(r, st%ndat)
end function force

! nrg_(1p|box)_* implementations that will use the potent_*, force_* previously defined 
#include "nrgimp.f90"

end module

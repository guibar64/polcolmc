module none
!! Null potential. Can be used for hard spheres.
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

  potent_r=0.D0
end function potent_r

real(8) function potent_r2(r2, nd)
  implicit none
  real(8), intent(in) :: r2
  type(NrgData), intent(in) :: nd

  potent_r2=0.D0
end function potent_r2


real(8) function force_r(r, nd)
  implicit none
  real(8), intent(in) :: r
  type(NrgData), intent(in) :: nd

  force_r=0.D0
end function force_r

real(8) function potential(r,i,j,st)
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: r
  type(PcmcState) :: st
  potential = potent_r(r, st%ndat)
end function potential

real(8) function force(r,i,j,st)
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: r
  type(PcmcState) :: st
  force = force_r(r, st%ndat)
end function force

! nrg_(1p|box)_* implementations that will use the potent_*, force_* previously defined 
#include "nrgimp.f90"

end module

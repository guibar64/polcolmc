module yukawa2
  !! Fast Yukawa potential
  use state

  private :: potent_r, potent_r2, force_r

! The nrg_(1p|box)_* should be private to the module
#include "nrgprivs.f90"

contains

subroutine initpot(st)
  use vscoul_mod
  use iso_fortran_env
  implicit none
  type(PcmcState) :: st
  real(8) :: r2min, r2max, rmin, rmax, kp
  logical :: vscoul_dbg
#ifdef DEBUG_YUK2
  vscoul_dbg=.true.
#else
  vscoul_dbg=.false.
#endif
  
  call gentab(r2min, r2max, vscoul_dbg)

  kp=1.D0/st%inp%debye_length
  st%ndat%kappa2 = kp*kp

  ! This method has its own separation range. Outside of it it breaks.
  rmin = sqrt(r2min)/st%ndat%inv_debye_length
  rmax = sqrt(r2max)/st%ndat%inv_debye_length
  if(rmin > minval(st%dist%rad)) then
     write(error_unit, '(A, G15.8, A, G15.8)') "yukawa2::initpot: Error: min separations (", rmin, &
          ") greater than min radius ", minval(st%dist%rad)
     stop 2
  end if
  ! For the upper bound, an early quit is not necessary, setting a maximum is enough.
  ! (This will be handled later in `init_modes_nrg_cutoff`.
  st%ndat%rcutoff_max = rmax

end subroutine initpot

real(8) function yuk2_core(p) result(t)
  use vscoul_mod
  real(8), intent(in) :: p
  integer(4) :: ixx(2)
  integer :: m
  real(8) :: xx
  equivalence (ixx, xx)
  xx = p
  !ixx = transfer(p,[0,0]) ! equivalence seems slighty faster

  m=iShft(ixx(imsw)-indx,-nShift)

  m=m*(nPol+1)
  t=c(m+4)
  t=t*p+c(m+3)
  t=t*p+c(m+2)
  t=t*p+c(m+1)
  t=t*p+c(m)
end function yuk2_core

real(8) function potent_r(r, nd)
  implicit none
  real(8), intent(in) :: r
  type(NrgData), intent(in) :: nd

  if(r<nd%rcutoff_max) then
    potent_r = nd%inv_debye_length*yuk2_core(r*r*nd%kappa2)
  else
    potent_r = exp(-r*nd%inv_debye_length)/r
  end if
end function potent_r


real(8) function potent_r2(r2, nd)
  implicit none
  real(8), intent(in) :: r2
  type(NrgData), intent(in) :: nd

  potent_r2 = nd%inv_debye_length*yuk2_core(r2*nd%kappa2)
  
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

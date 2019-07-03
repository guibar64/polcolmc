#define FFTW3 3
#define SLOWFFT 0

#ifndef FOURIER_TRANSFORM
#define FOURIER_TRANSFORM SLOWFFT
#endif

module sinft
!! Discrete sinus transform
!!
!! By default, fftw3 routines (not included) are not used (this can be changed by the appropriate flags
!! at compilation).
!! They are way faster, but for now FT is only used at the end of a simulation, the speed-up
!! is negligible for a real simulation.
use iso_c_binding
implicit none

real(8), parameter, private :: PI=3.141592653589793

! interface to fftw3
integer(C_INT), parameter :: FFTW_RODFT10 = 9,FFTW_ESTIMATE = 64
interface
  function fftw_plan_r2r_1d(n,in,out,kind,flags) bind(C)
    import
    type(c_ptr) fftw_plan_r2r_1d
    integer(c_int), value :: n, flags
    real(c_double) :: in(*),out(*)
    integer(4), value :: kind
  end function

  subroutine fftw_execute_r2r(p,in,out) bind(c)
    import
    type(c_ptr), value :: p
    real(c_double) :: in(*),out(*)
  end subroutine

  subroutine fftw_destroy_plan(p) bind(c)
    import
    type(c_ptr), value :: p
  end subroutine
end interface

contains

#if FOURIER_TRANSFORM == FFTW3
subroutine sinft_calc(n,in,out)
  integer, intent(in) :: n
  real(8) :: in(n), out(n)
  if(nfft /= n) then
    ! Plan is reused as long as data length is the same. ⚠ Not thread safe.
    if(nfft /= 0)   call fftw_destroy_plan(pfft)
    pfft = fftw_plan_r2r_1d(n, in, out, FFTW_RODFT10, FFTW_ESTIMATE)
  end if
  call fftw_execute_r2r(pfft,in,out)
end subroutine

subroutine sinft_cleanup()
  if(nfft /= 0)   call fftw_destroy_plan(pfft)
end subroutine

#else
subroutine sinft_calc(n, in, out)
  !! Computes a discrete sinus transform.
  !! @Warning
  !! input dta `in` is mutated.
  !! @endwarning
  integer, intent(in) :: n !! data length
  real(8),intent(inout) :: in(n) !! input data
  real(8), intent(out) :: out(n) !! result

  integer :: k,j
  real(8) :: sj,o,pion
  complex(8) :: w00, w0, ww
  pion=PI/real(n,8)
  w00=exp(cmplx(0.d0,pion,kind=8))
  w0=cmplx(1.d0,0.d0,kind=8)
  ww=w0
  do k=0,n-1
    o = 0.d0
    w0 = w0 * w00
    ww = sqrt(w0)
    do j=0,n-1
      sj = aimag(ww)
      o = o + in(j+1)*sj
      ww = ww * w0
    end do
    out(k+1) = o
  end do
  out(:) = out(:) * 2
end subroutine

subroutine sinft_cleanup()
!! Cleanup fft "plans". Only relevant for fftw3.
end subroutine
#endif

end module
module struct_fact
!! Computes the struture factor from the radial distribution functions.
use iso_fortran_env
use iso_c_binding
use sinft
implicit none

real(8), parameter, private :: PI=3.141592653589793

private :: eleFFS, do_dst

contains

function EleFFS(q, r)
  !! Proportional to the form factor of a sphere.
  real(8) :: EleFFS
  real(8), intent(in) :: q,r
  real(8) :: x
  x=q*r
  eleFFs = sin(x)-x*cos(x)
end function

subroutine do_dst(n, y, out)
  !! Do a sinus tranform of y-1.
  integer :: n !! dat length
  real(8), intent(in) :: y(n) !! input data
  real(8) :: out(n) !! result
  real(8) :: in(n)
  integer :: i
  do i=1,n
    in(i) = (y(i)-1.D0)*(real(i,8)-0.5D0)
  end do
  call sinft_calc(n, in, out)
end subroutine

subroutine calc_sdeq( dx, n,  q,  sf,  pf, nqext,qext,  sfext, ng,  gder, nf,  pop&
	      &,  rad,  phi)
  !! Calculates S(q) from RDFs.
  use spline
  integer, intent(in) :: n !! length of q array
  integer, intent(in) :: nqext !! length of qext
  integer, intent(in) :: ng  !! lenght of RDFs
  integer, intent(in) :: nf  !! number of families
  real(8), intent(in) :: dx  !! space bin of RDFs
  real(8), intent(in) :: phi !! volume fraction
  real(8) :: q(n) !! wave numbers
  real(8) :: sf(n) !! struture factor
  real(8) :: pf(n) !! form factor
  real(8) :: qext(nqext) !! new q values
  real(8) :: sfext(nqext) !! interpolated struture factor
  real(8) :: gder(n,ng)  !! RDFs
  real(8) :: pop(nf)     !! family populations
  real(8) :: rad(nf)     !! family radii
  
  integer :: i, j,k,ind
  real(8) :: vp, gsk, dump,dsq
  real(8), allocatable :: s(:),g(:),rho(:),sf2ip(:)
  allocate(s(n),g(n),rho(nf))

  if(ng /= (nf*(nf+1))/2) stop "calc_sqdeq: ng != (nf*(nf+1))/2"

  vp=0.0
  do i=1,nf
    vp = vp + pop(i)*4.D0/3.D0*PI*(rad(i)**3)
  end do
  do i=1,nf
    rho(i) = pop(i)*phi/vp
  end do
  
  do k=1,n
    q(k) = (PI*(k))/(dx*n);
  end do

  sf = 0.D0
  pf = 0.D0
  
  ind = 0
  do i=1,nf
    ind = ind + 1
    do k=1,n
      g(k) = gder(k,ind)
    end do
    call do_dst(n, g, s)
    do k=1,n
      gsk = 1.d0 + (2.d0*PI*(rho(i)*dx*dx)/q(k))*s(k)
      dump = EleFFS(q(k), rad(i))*EleFFS(q(k), rad(i))*rho(i)
      sf(k) = sf(k) + gsk*dump
      pf(k) = pf(k) + dump
    end do
  end do
  do i=1,nf
    do j=i+1,nf
      ind = ind + 1
      do k=1,n
        g(k) = gder(k,ind)
      end do
      call do_dst(n, g, s)
      do k=1,n
        dsq = sqrt(rho(i)*rho(j))
        gsk = 2.D0*PI*(dsq*dx*dx)/q(k)*s(k)
        dump = EleFFS(q(k), rad(i))*EleFFS(q(k), rad(j))*dsq;
        sf(k) = sf(k) + 2.D0*gsk*dump
      end do
    end do
  end do
  call sinft_cleanup()
  sf = sf / pf

  allocate(sf2ip(n))
  call cubspline_init(n, q, sf,sf2ip)
  call cubspline_equalsep_eval(n,q,sf,sf2ip,nqext,qext, sfext)

end subroutine

end module struct_fact

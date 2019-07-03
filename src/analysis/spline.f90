module spline
!! Spline utilities
implicit none

private :: findIdx

contains

subroutine cubspline_init(n, x, y,y2)
  !! Initialize a cubic spline
  integer, intent(in) :: n !! data length
  real(8), intent(in) :: x(n) !! x input data
  real(8), intent(in) :: y(n) !! y(x) 
  real(8), intent(out) :: y2(n) !! spline data

  integer :: i,k
  real(8) :: p,qn,sig, un,u(n), yp1,yp2,h

  h = x(2)-x(1)

  y2(1)=0.d0
  u(1)=0.d0
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    yp1=(y(i+1)-y(i))/(x(i+1)-x(i))
    yp2=(y(i)-y(i-1))/(x(i)-x(i-1))
    u(i)=(6.d0*(yp1-yp2)/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  qn=0.d0
  un=0.d0
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
  do k=n-1,1,-1
    y2(k) = y2(k)*y2(k+1) + u(k)
  end do
  
end subroutine

function cubspline_equalsep_eval_single_nchk(hi, h2sur6,xmin,y,y2, val) result(res)
  !! Evaluate a spline of equally sperated absissa at `val`. Input range is not checked. Helper function.
  real(8) :: res
  real(8), intent(in) :: y(*) !! y(x) data points
  real(8), intent(in) :: y2(*) !! spline data
  real(8), intent(in) :: val  !! value
  real(8), intent(in) :: hi !! 1/h,  where h is step
  real(8), intent(in) :: h2sur6 !! h^2/6, where h is step
  real(8), intent(in) :: xmin !! minimum x value
  integer :: nr2
  real(8) :: a,b

  a=(val-xmin)*hi
  nr2=int(a)
  a = a - real(nr2,8)
  b=1.d0-a
  res=y(nr2+2)*a+y(nr2+1)*b+&
    &(y2(nr2+2)*a*(a*a-1.d0)+y2(nr2+1)*b*(b*b-1.d0))*h2sur6
end function

subroutine cubspline_equalsep_eval(n, x, y, y2, nnew,xnew, ynew)
  !! Evaluates an array of values.
  integer, intent(in) :: n !! data length
  real(8), intent(in) :: x(n) !! input values
  real(8), intent(in) :: y(n) !! y(x) data values
  real(8), intent(in) :: y2(n) !! spline data
  integer, intent(in) :: nnew !! output data length
  real(8), intent(in) :: xnew(nnew) !! values to be evaluated
  real(8), intent(out) :: ynew(nnew) !! result array y(xnew)
 
  integer :: i,imin, imax
  real(8) :: h2s6,hi,h
  h = x(2) - x(1)
  h2s6 = h*h/6.d0
  hi = 1.d0 / h

  do imin=1,nnew
    if(xnew(imin) > x(1)) exit
    ynew(imin) = y(1)
  end do
  do imax=nnew,1,-1
    if(xnew(imax) < x(n)) exit
    ynew(imax) = y(n)
  end do

  do i=imin,imax
    ynew(i) = cubspline_equalsep_eval_single_nchk(hi, h2s6,x(1),y,y2, xnew(i))
  end do 

end subroutine


integer function findIdx(n, x, val)
  !todo: store klo, khi for sequential calls
  integer, intent(in) :: n
  real(8), intent(in) :: x(n)
  real(8), intent(in) :: val
  integer :: klo, khi,k
  klo = 1
  khi = n
  do while(khi-klo>1)
    k=(khi+klo) / 2
    if(x(k) > val) then
      khi = k
    else
      klo=k
    end if
  end do
  findIdx= klo
end function

real(8) function cubspline_eval(n, x, y, y2, val)
  !! Evaluates a value. `val` must be in range.
  integer, intent(in) :: n !! data length
  real(8), intent(in) :: x(n) !! input values
  real(8), intent(in) :: y(n) !! y(x) data values
  real(8), intent(in) :: y2(n) !! spline data
  real(8), intent(in) :: val !! value
 
  integer :: nr2
  real(8) :: a,dx,b,c,d,h2sur6

  nr2 = findIdx(n,x, val)
  dx = x(nr2+1) - x(nr2)
  a = (val - x(nr2))/dx
  b = 1.0d0 - a
  c = a*(a*a-1.0d0)
  d = b*(b*b - 1.0d0)
  h2sur6 = dx*dx*(1.0d0/6.0d0)
  cubspline_eval = y(nr2 + 1) * a + y(nr2) * b + (y2(nr2 + 1) * c + y2(nr2) * d )* h2sur6


end function

end module

module SpherHarm
  !! Functions related to spherical harmonics
  implicit none
contains
  
complex(8) function SpherHarmSpecial(l,m,un)
  !! Ylm as a function of position vector
  integer,intent(in) :: l, m
  real(8),intent(in) :: un(3)
  real(8), parameter :: PI=3.14159265358979D0
  integer :: i,fact,mm
  real(8) ::  phi
  mm=abs(m)
  
  fact=1
  do i=l+mm,l-mm+1,-1
     fact=fact*i
  end do
  
  phi=atan2(un(2),un(1))

  SpherHarmSpecial=sqrt((2*l+1)/(4*fact*PI))*plgndr(l,mm,un(3))*&
      cmplx(cos(m*phi),sin(m*phi), 8)
end function SpherHarmSpecial


real(8) function SpherHarmTrunc(l,m,costheta)
  !! Ylm(θ,ϕ) without the exp(i m ϕ) factor 
  integer,intent(in) :: l,m
  real(8),intent(in) :: costheta  ! cos(θ)
  real(8), parameter :: PI=3.14159265358979D0
  integer :: i,fact,mm
  mm=abs(m)
  
  fact=1
  do i=l+mm,l-mm+1,-1
     fact=fact*i
  end do

  SpherHarmTrunc=sqrt((2*l+1)/(4*fact*PI))*plgndr(l,mm,costheta)
end function SpherHarmTrunc

function plgndr(l,m,x)
  !! Computes the associated renormalized Legendre polynomial Pm
  !! l (x). m and l must statisfy 0<=m<=l and x must be in [-1,1].
  !! Legendre polynomial.
  integer, intent(in) :: l,m
  real(8), intent(in) :: x
  real(8) :: plgndr 
  integer :: ll
  real(8) :: fact,pll,pmm,pmmp1,somx2
  pll=0.d0
  pmm=1.0 !Compute Pmm
  if (m > 0) then
     somx2=sqrt((1.0D0-x)*(1.0D0+x))
     fact=1.
     do ll=1,m
        pmm=-pmm*fact*somx2
        fact=fact+2.
     enddo
  endif
  
  if (l == m) then
     plgndr=pmm
  else
     pmmp1=x*(2*m+1)*pmm !Compute Pmm+1.
     if (l == m+1) then
        plgndr=pmmp1
     else !Compute Pml , l > m+ 1.
        do ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm=pmmp1
           pmmp1=pll
        end do
        plgndr=pll
     end if
  end if
end function plgndr

end module

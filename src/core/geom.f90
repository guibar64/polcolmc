module geom
  !! Some operations on box tensors etc...
  implicit none

contains

  real(8) function inverse_box(a,ai)
    !! Computes the inverse of box `a` to `ai`. Returns the determinant of `a`.
    implicit none
    real(8) :: a(3,3),ai(3,3)

    real(8) :: det

    det = +a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))&
      -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))&
      +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))


    if(det>0.d0) then
      ai(1,1)=+(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      ai(1,2)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      ai(1,3)=+(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      ai(2,1)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      ai(2,2)=+(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      ai(2,3)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      ai(3,1)=+(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      ai(3,2)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))
      ai(3,3)=+(a(1,1)*a(2,2)-a(1,2)*a(2,1))
      
      ai = ai/det
    else
      ai = 0.d0
    end if

    inverse_box = det

  end function inverse_box

  subroutine print_box(unit, a)
    !! Pints a 3x3 array `a` to `unit` in a line1\nline2\n... fashion.
    implicit none
    integer :: unit
    real(8) :: a(3,3)
    integer :: i
    do i=1,3
      write(unit, '(3G16.7)') a(i,1:3)
    end do
  end subroutine print_box

  real(8) function volume_box(a)
    !! Compute the volume of box of tensor `a`.
    implicit none
    real(8) :: a(3,3)

    real(8) :: det


    det=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) - a(2,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3))&
      +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))

    volume_box = det

  end function volume_box

  subroutine calc_metric(a, m)
    !! Compute the metric of a box of tensor `a` to `m`.
    real(8), intent(in) :: a(3,3)
    real(8), intent(out) :: m(3,3)

    integer :: i,j,k
    do i=1,3
      do j=1,3
        m(i,j) = 0.d0
        do k=1,3
          m(i,j) = m(i,j) + a(k,i)*a(k,j)
        end do
      end do
    end do

  end subroutine calc_metric

  function matvec3_prod(bt,p)
    !! Matrix product of `bt` with `p`.
    real(8) :: matvec3_prod(3)
    real(8), intent(in) :: bt(3,3), p(3)
    integer :: i
    do i=1,3
      matvec3_prod(i) = bt(i,1)*p(1) + bt(i,2)*p(2) + bt(i,3)*p(3)
    end do
  end function matvec3_prod

  subroutine create_box_from_3l3a(boxtens,La,Lb,Lc,alpha,beta,gamma)
    !! Computes a box tensor from box lengths and box angles.
    implicit none
    real(8),intent(in) :: La,Lb,Lc,alpha,beta,gamma
    real(8),intent(out) :: boxtens(3,3)

    real(8) :: a(3),b(3),c(3)

    a(1)=la ; a(2)=0.D0 ; a(3)=0.D0
    b(1)=lb*cos(gamma) ; b(2) = lb*sin(gamma) ; b(3)=0.D0
    c(1)=lc*cos(beta) 
    c(2)=lc/(lb*sin(gamma))*(lb*cos(alpha)-la*cos(beta)*cos(gamma))
    c(3) = sqrt(lc**2-c(2)**2-c(1)**2)

    boxtens(1,:) = a
    boxtens(2,:) = b
    boxtens(3,:) = c

  end subroutine create_box_from_3l3a

end module geom
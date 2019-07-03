elemental function pbc(x) result(res)
  !! returns coordinate difference of first minimum image
  real(8)  :: res
  real(8), intent(in) :: x
  res = x - int(2._8*x)
end function pbc

real(8) function dist2_mi_generic(pos,posout, MetricG)
  !! Returns square of the separation between the minimum image of two particles
  !! for a generic (triclinic) box.
  !! `posout` is set to the vector between the images.
  use iso_fortran_env
  implicit none
  real(8),intent(in) :: pos(3)
  real(8),intent(out) :: posout(3)
  real(8),intent(in) :: MetricG(3,3)
  
  integer :: j

  do j=1,3
     posout(j) = pbc(pos(j))
  end do
  
  dist2_mi_generic = posout(1)*(posout(1)*MetricG(1,1) + 2*posout(2)*MetricG(2,1)+2*posout(3)*MetricG(3,1))+&
       posout(2)*(posout(2)*MetricG(2,2) + 2*posout(3)*MetricG(3,2))+&
       posout(3)*(posout(3)*MetricG(3,3))
  
end function dist2_mi_generic

real(8) function dist2_mi_cubic(pos,posout,MetricG)
  !! Returns square of the separation between the minimum image of two particles
  !! for a cubic box.
  !! `posout` is set to the vector between the images.
  implicit none
  real(8),intent(in) :: pos(3)
  real(8),intent(out) :: posout(3)
  real(8),intent(in) :: MetricG(3,3)
  
  integer :: j

  do j=1,3
     posout(j) = pbc(pos(j))
  end do
  dist2_mi_cubic=posout(1)*posout(1)+posout(2)*posout(2)+posout(3)*posout(3)
  dist2_mi_cubic = dist2_mi_cubic*MetricG(1,1)
  
end function dist2_mi_cubic

real(8) function dist2_mi_cubicXY(pos,posout,MetricG)
  !! Returns square of the separation between the minimum image of two particles
  !! for a cubic box with periodic boundary conditions only in the XY directions.
  !! `posout` is set to the vector between the images.
  implicit none
  real(8),intent(in) :: pos(3)
  real(8),intent(out) :: posout(3)
  real(8),intent(in) :: MetricG(3,3)
  integer :: j

  do j=1,2
     posout(j) = pbc(pos(j))
  end do

  posout(3)=pos(3)

  dist2_mi_cubicXY=posout(1)*posout(1)+posout(2)*posout(2)+posout(3)*posout(3)
  dist2_mi_cubicXY = dist2_mi_cubicXY*MetricG(1,1) 
 
  
end function dist2_mi_cubicXY

real(8) function dist2_mi_hexagonal(pos,posout, MetricG)
  !! Returns square of the separation between the minimum image of two particles
  !! for an hexagonal box.
  !! `posout` is set to the vector between the images.
  use iso_fortran_env
  implicit none
  real(8),intent(in) :: pos(3)
  real(8),intent(out) :: posout(3)
  real(8),intent(in) :: MetricG(3,3)
  
  integer :: j

  do j=1,3
     posout(j) = pbc(pos(j))
  end do
  dist2_mi_hexagonal = MetricG(1,1)*(posout(1)*posout(1) + posout(2)*(posout(2) + posout(1)))+&
       posout(3)*posout(3)*MetricG(3,3)
  
end function dist2_mi_hexagonal


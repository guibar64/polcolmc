module histograms
  implicit none
  !! Simple histogram types and methods
  
  type HistoR8
    integer(8), allocatable :: nh(:)
    real(8), allocatable :: h(:)
    real(8) :: dri, rmin
    integer(8) :: nupd
  end type

contains

subroutine histor8_init(p, dr, rmin, n)
  !! Initialize an histogram
  type(HistoR8), intent(out) :: p
  real(8), intent(in) :: dr   !! step
  real(8), intent(in) :: rmin !! minimum value
  integer, intent(in) :: n                !! number of points 

  p%nupd = 0
  p%dri = 1.d0 / dr
  p%rmin = rmin
  if(allocated(p%nh)) then
    if(ubound(p%nh, 1) /= n) then
      deallocate(p%nh)
      allocate(p%nh(n))
    end if
  else
    allocate(p%nh(n))
  end if
  p%nh = 0.0

end subroutine

subroutine histor8_reset(p)
  !! Reset an histogram to do a new sampling
  type(HistoR8), intent(inout) :: p

  p%nupd = 0
  p%nh = 0.0
end subroutine

subroutine histor8_update(p, val)
  !! Update an histogram
  type(HistoR8), intent(inout) :: p
  real(8) :: val
  integer :: ival
  ival = floor(p%dri * (val-p%rmin)) + 1
  if(ival >= lbound(p%nh,1) .and. ival <= ubound(p%nh,1)) then
    p%nh(ival) = p%nh(ival) + 1
  end if
  p%nupd = p%nupd + 1

end subroutine

subroutine histor8_update_n(p, val, n)
  !! Update an histogram with `n` values
  type(HistoR8), intent(inout) :: p
  real(8), intent(in) :: val
  integer, intent(in) :: n
  integer :: ival
  ival = floor(p%dri * (val-p%rmin)) + 1
  if(ival >= lbound(p%nh,1) .and. ival <= ubound(p%nh,1)) then
    p%nh(ival) = p%nh(ival) + n
  end if
  p%nupd = p%nupd + n
end subroutine


subroutine histor8_set_distribution(p)
  !! Sets p%h to a distribution
  type(HistoR8), intent(inout) :: p

  real(8) :: toti

  if(allocated(p%h)) then
    if(ubound(p%h,1) /= ubound(p%nh,1)) then
      deallocate(p%h)
      allocate(p%h(ubound(p%nh,1)))
    end if
  else
    allocate(p%h(ubound(p%nh,1)))
  end if

  toti = 1.d0/sum(p%nh) ! or p%nupd 

  p%h = real(p%nh,8)*p%dri*toti

end subroutine
    
end module
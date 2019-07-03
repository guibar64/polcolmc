
function nrg_1p(st,b,par)
	!! Computes energy of one particle
  type(NrgRec) :: nrg_1p
  type(PcmcState), intent(in) :: st
  type(Box),intent(in),target :: b   !! simulation box
  type(Particule), intent(in) :: par !! the particle
  select case(st%mode_combo)
{0}
    case default
      stop 9
  end select
end function

function nrg_box(st,b,press)
	!! Computes energy of a box.
	!! If `press` is called, computes the pressure and returns it in `press`.
  type(NrgRec) :: nrg_box
  type(PcmcState), intent(in) :: st
  type(Box), intent(in) :: b                     !! simulation box
  type(Pressure),intent(out), optional :: press  !! Pressure
  select case(st%mode_combo)
{1}
    case default
      stop 9
  end select
end function nrg_box

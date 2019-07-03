include "dist2_mi.f90"


function nrg_1p(st,b,par)
  !! Computes energy of one particle
  type(NrgRec) :: nrg_1p
  type(PcmcState), intent(in) :: st
  type(Box),intent(in),target :: b   !! simulation box
  type(Particule), intent(in) :: par !! the particle
  select case(st%mode_combo)
    case(0)
      nrg_1p = nrg_1p_standard_generic(st,b,par)
    case(1)
      nrg_1p = nrg_1p_standard_cubic(st,b,par)
    case(2)
      nrg_1p = nrg_1p_standard_cubicXY(st,b,par)
    case(3)
      nrg_1p = nrg_1p_standard_hexagonal(st,b,par)
    case(100)
      nrg_1p = nrg_1p_celldec_generic(st,b,par)
    case(101)
      nrg_1p = nrg_1p_celldec_cubic(st,b,par)
    case(102)
      nrg_1p = nrg_1p_celldec_cubicXY(st,b,par)
    case(103)
      nrg_1p = nrg_1p_celldec_hexagonal(st,b,par)

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
    case(0)
      nrg_box = nrg_box_standard_generic(st,b,press)
    case(1)
      nrg_box = nrg_box_standard_cubic(st,b,press)
    case(2)
      nrg_box = nrg_box_standard_cubicXY(st,b,press)
    case(3)
      nrg_box = nrg_box_standard_hexagonal(st,b,press)
    case(100)
      nrg_box = nrg_box_celldec_generic(st,b,press)
    case(101)
      nrg_box = nrg_box_celldec_cubic(st,b,press)
    case(102)
      nrg_box = nrg_box_celldec_cubicXY(st,b,press)
    case(103)
      nrg_box = nrg_box_celldec_hexagonal(st,b,press)

    case default
      stop 9
  end select
end function nrg_box

#define DIST2_MI dist2_mi_generic
#define NRG_1P  nrg_1p_standard_generic
#define NRG_BOX  nrg_box_standard_generic
#include "nrg_standard.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_cubic
#define NRG_1P  nrg_1p_standard_cubic
#define NRG_BOX  nrg_box_standard_cubic
#include "nrg_standard.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_cubicXY
#define NRG_1P  nrg_1p_standard_cubicXY
#define NRG_BOX  nrg_box_standard_cubicXY
#include "nrg_standard.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_hexagonal
#define NRG_1P  nrg_1p_standard_hexagonal
#define NRG_BOX  nrg_box_standard_hexagonal
#include "nrg_standard.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_generic
#define NRG_1P  nrg_1p_celldec_generic
#define NRG_BOX  nrg_box_celldec_generic
#include "nrg_celldec.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_cubic
#define NRG_1P  nrg_1p_celldec_cubic
#define NRG_BOX  nrg_box_celldec_cubic
#include "nrg_celldec.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_cubicXY
#define NRG_1P  nrg_1p_celldec_cubicXY
#define NRG_BOX  nrg_box_celldec_cubicXY
#include "nrg_celldec.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

#define DIST2_MI dist2_mi_hexagonal
#define NRG_1P  nrg_1p_celldec_hexagonal
#define NRG_BOX  nrg_box_celldec_hexagonal
#include "nrg_celldec.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX

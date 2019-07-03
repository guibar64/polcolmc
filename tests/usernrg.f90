
program test_usernrg
  use simulation
  implicit none
  type(MCSimulation) :: simu

  ! Always before call to 'init_simulation'
  call set_nrg_routines(simu%nrg, mypotent, myforce, mynrg_1p, mynrg_box, myinitpot)
  call init_simulation(simu, "distrib", "usernrg.cfg", "usernrg_in.inst")
  call run_simulation_equil(simu, 20)
  call run_simulation(simu, 10)

  call check_double(simu%state%res(1)%ep, 26221.73d0)
  call check_double(simu%smax(1), 2.320693d0)
  
contains

  subroutine myinitpot(st)
    type(PcmcState) :: st
  end subroutine myinitpot

  real(8) function myforce(r,i,j,st)
    integer, intent(in) :: i, j
    real(8), intent(in) :: r
    type(PcmcState) :: st
    myforce = 0.d0
  end function myforce

  real(8) function mypotent(r,i,j,st)
    integer, intent(in) :: i, j
    real(8), intent(in) :: r
    type(PcmcState) :: st
    mypotent = 0.d0
  end function mypotent

  real(8) function corepot(pos, met)
    real(8) :: pos(3), met(3,3)
    real(8) :: relpos(3)
    relpos = pos - [0.5d0, 0.5d0, 0.5d0]
    corepot = dot_product(relpos, matvec3_prod(met, relpos)) * 1.d-4
  end function corepot
  
  function  mynrg_1p(st,b,par)
    type(NrgRec) :: ener, mynrg_1p
    type(PcmcState), intent(in) :: st
    type(Box),intent(in),target :: b
    type(Particule), intent(in) :: par

    ! basically check overlaps
    ener = none_nrg_1p(st,b,par)
    if(ener%noverlap) then
      mynrg_1p%noverlap = .true.
      mynrg_1p%ep = corepot(par%pos, b%met)
    else
      mynrg_1p%noverlap = .false.
    end if
    
  end function mynrg_1p

  function mynrg_box(st,b,press)
    type(NrgRec) :: mynrg_box
    type(PcmcState), intent(in) :: st
    type(Box), intent(in) :: b
    type(Pressure),intent(out),optional :: press

    type(NrgRec) :: te
    logical :: over
    real(8) :: e
    integer :: i
    ! basically check overlaps
    te = none_nrg_box(st,b,press)
    if(te%noverlap) then
      mynrg_box%noverlap = .true.
      e = 0.d0
      do i=1,b%N
        e = e + corepot(b%parts(i)%pos, b%met)
      end do
      mynrg_box%ep = e
    else
      mynrg_box%noverlap = .false.
    end if

  end function mynrg_box
  
  include "check.f90"

end program test_usernrg


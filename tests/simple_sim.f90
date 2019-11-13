
program simple_sim
  use simulation
  use iso_fortran_env
  implicit none
  type(MCSimulation) :: simu
  real(8), parameter :: ep=7747.465d0,press=1947.6201d0
  integer :: step
  call init_simple_simulation(simu,"distrib","simple.cfg","simple_in.inst")
  call init_step_series(simu, "pcmc_simple_traj.inst","pcmc_simple_")
  do step=1,30
      call step_series_next(simu, step)
  end do
  call step_series_end(simu)

  call check_double(simu%state%res(1)%ep, ep)
  call check_double(simu%state%res(1)%P, press )

  contains

  include "check.f90"

end program

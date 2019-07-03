
program gibbst27
  use simulation
  implicit none
  type(mcsimulation) :: simu
  real(8) :: press(2) = [22404.733,22379.507]
  real(8) :: ep(2) = [2572.343, 2587.491]
  real(8) :: vf(2) = [0.21897455, 0.22098973]

  call init_simulation(simu,"dist-gibbst27","gibbst27.cfg","gibbst27_in.inst")
  call run_simulation(simu, 100)
  call check_double(simu%state%res(1)%P, press(1))
  call check_double(simu%state%res(1)%Ep, ep(1))
  call check_double(simu%state%res(1)%hetv, vf(1))

  call check_double(simu%state%res(2)%P, press(2))
  call check_double(simu%state%res(2)%Ep, ep(2))
  call check_double(simu%state%res(2)%hetv, vf(2))

  contains

  include "check.f90"

end program

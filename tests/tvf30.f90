
program tfv30
  use simulation
  implicit none
  call test1()
  call test2()

contains

  include "check.f90"

  subroutine test1()
    type(mcsimulation) :: simu

    call init_simulation(simu,"distrib","vf30.cfg","in.inst")
    call run_simulation_equil(simu, 10)
    call run_simulation(simu, 10)

    call check_double(simu%state%res(1)%P, 9485.387d0)
    call check_double(simu%state%res(1)%Ep, 31546.10d0)
    call check_double(simu%Smax(1), 1.6220679d0)
  end subroutine test1

  subroutine test2()
    type(mcsimulation) :: simu

    call init_simulation(simu,"distrib","vf30_2.cfg","in.inst")
    call run_simulation_equil(simu, 10)
    call run_simulation(simu, 10)

    call check_double(simu%state%res(1)%P, 9432.0292d0)
    call check_double(simu%state%res(1)%Ep, 31804.09d0)
    call check_double(simu%Smax(1), 1.6201629d0)
  end subroutine test2

end program tfv30


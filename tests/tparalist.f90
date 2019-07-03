
program simple_sim
  use simulation
  use readparam
  use iso_fortran_env
  implicit none
  integer :: i
  real(8) :: rc
  logical :: prod
  type(MCSimulation) :: simu
  character(len=:), allocatable :: val

  call init_simple_simulation(simu,"distrib","vf30.cfg","simple_in.inst")
  
  if(get_param(simu%state%config, "log_period", val)) then
    read(val,*) i
    if(i /= 10000) then
      write(*,*) "log_period /= ", 10000
      stop 1
    endif
  else
    write(*,*) "'log_period' not found"
    stop 1
  end if

  if(get_param(simu%state%config, "equilibration_steps", val)) then
    read(val,*) i
    if(i /= 10) then
      write(*,*) "equilibration_steps /= ", 10000
      stop 1
    endif
  else
    write(*,*) "'equilibration_steps' not found"
    stop 1
  end if

  if(get_param(simu%state%config, "production", val)) then
    prod = read_logical(val)
    if(.not. prod) then
      write(*,*) "production /= ", .true.
      stop 1
    endif
  else
    write(*,*) "'production' not found"
    stop 1
  end if

  if(get_param(simu%state%config, "bondorder.cutoff", val)) then
    read(val, *) rc
    if(rc /= 55.8_8) then
      write(*,*) "key /= ", 55.8
      stop 1
    endif
  else
    write(*,*) "'bondorder.cutoff' not found"
    stop 1
  end if

  call warn_about_unkown_params(output_unit, simu%state%config)

  contains

  include "check.f90"

end program

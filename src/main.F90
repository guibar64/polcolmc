
program pdcolmc
!! main program
use simulation, only: MCSimulation, run_simulation_equil, run_simulation, init_simulation, set_threads_openmp, &
  deallocate_simulation
use tools
use readparam, only: warn_about_unkown_params
use sgc, only: sgc_doublerun
use iso_fortran_env
!$ use omp_lib
implicit none

character(*), parameter :: POLCOLMC_VERSION="0.7.0"

integer, parameter :: RUN_SIMU=0, CALC_GDER=1, CONVERT=2, RESCALE=3, EXTRACT=4, DOUBLE=5,TABLES=6, CALC_BONDORDER=7,&
  CALC_MSD=8, PRUNE=9

logical :: para_lec_only=.false.,verb=.false., use_bocfg=.false., bo_multiconf = .false.
character(256) :: option
character(256) :: progname,fich_dist="distrib",fich_par="polcolmc.cfg",fich_in="in.inst",fich_out="polcolmc.out",&
  fich_bocfg = "bondorder.cfg"
character(256) :: infile, outfile,iformat="",oformat=""
character(256) :: infiles(43)
integer :: ni, partselected=1, action=RUN_SIMU, ibox=1, ninfiles=0, bstep=0, estep=huge(1), period_pr = 1
integer, allocatable :: ssseed(:)

integer :: i,argc
real(8) :: dump, factor=1.d0, radius_extract = 10.0,dr_crdf=0.1,dr_mktables=0.1D0

integer :: nthreads_req=0

type(MCSimulation) :: the_sim

! Parse comand line.
call get_command_argument(0, PROGNAME)
argc=command_argument_count()
i=0
do while (i < argc)
  i=i+1
  call get_command_argument(i,option)
  select case(trim(option))
  case("-v")
    verb=.true.
  case("--version")
    call write_version()
    stop 0
  case("-l")
    para_lec_only=.true.
  case("-i")
    i=i+1
    call get_command_argument(i,fich_in)
  case("-d")
    i=i+1
    call get_command_argument(i,fich_dist)
  case("-p")
    i=i+1
    call get_command_argument(i,fich_par)
  case("mktables")
    action = TABLES
    i=i+1 ; call check_arg(i, "BIN")
    call get_command_argument(i,option) 
    read(option,*) dr_mktables
  case("--nthreads")
    i=i+1
    call get_command_argument(i,option)
    read(option,*) nthreads_req
  case("-h","--help")
    call write_usage()
    stop 0
  case("run")
    action = RUN_SIMU
  case("rdf")
    action = CALC_GDER
    i=i+1 ; call check_arg(i, "INPUT_FILE")
    call get_command_argument(i,infile)
    i=i+1  ; call check_arg(i, "OUTPUT_FILE")
    call get_command_argument(i,outfile)
  case("--bstep")
    i=i+1 ; call check_arg(i, "BSTEP")
    call get_command_argument(i,option)
    read(option,*) bstep
  case("--estep")
    i=i+1 ; call check_arg(i, "ESTEP")
    call get_command_argument(i,option)
    read(option,*) estep
  case("rescale")
    action = RESCALE
    i=i+1 ; call check_arg(i, "INPUT_FILE")
    call get_command_argument(i,infile)
    i=i+1 ; call check_arg(i, "OUTPUT_FILE")
    call get_command_argument(i,outfile)
    i=i+1 ; call check_arg(i, "FACTOR")
    call get_command_argument(i,option)
    read(option,*) factor
  case("extract")
    action = EXTRACT
    i=i+1 ; call check_arg(i, "INPUT_FILE")
    i=i+1 ; call check_arg(i, "OUTPUT_FILE")
    call get_command_argument(i,outfile)
  case("--radius")
    i=i+1  ; call check_arg(i, "RADIUS")
    call get_command_argument(i,option)
    read(option,*) radius_extract
  case("--particle")
    i=i+1 ; call check_arg(i, "PARTICLE")
    call get_command_argument(i,option)
    read(option,*) partselected
  case("--dr")
    i=i+1  ; call check_arg(i, "DR")
    call get_command_argument(i,option)
    read(option,*) dr_crdf
  case("double")
    action = DOUBLE
    i=i+1 ; call check_arg(i, "INPUT_FILE")
    call get_command_argument(i,infile)
    i=i+1 ; call check_arg(i, "OUTPUT_FILE")
    call get_command_argument(i,outfile)
  case("convert")
    action = CONVERT
    
    ninfiles = 0
    do while (i < argc)
      i=i+1 ; call get_command_argument(i,infile)
      if(infile(1:1) == '-') then ! option, handled by main select
        i = i - 1 ! Whys is it needed
        exit
      end if
      ninfiles = ninfiles + 1
      infiles(ninfiles) = infile
    end do
    if(ninfiles < 2) then
      write(*,'(A)') "Syntax Error: needs at least two arguments: INPUT_FILE OUTPUT_FILE"
      stop 1
    end if
  case("--if")
    i = i + 1 ; call check_arg(i, "FORMAT")
    call get_command_argument(i, iformat)
  case("--of")
    i = i + 1 ; call check_arg(i, "FORMAT")
    call get_command_argument(i, oformat)
  case("--box")
    i=i+1  ; call check_arg(i, "BOX")
    call get_command_argument(i,option)
    read(option,*) ibox
  case("--bondorder-cfg", "--bo-cfg")
    i=i+1 ; call check_arg(i, "BOND_ORDER_CONFIG")
    call get_command_argument(i,option)
    fich_bocfg = trim(option)
    use_bocfg = .true.
  case("--multiconf")
    bo_multiconf = .true.
  case("bondorder")
    action = CALC_BONDORDER
    i=i+1  ; call check_arg(i, "INPUT_FILE")
    call get_command_argument(i,infile)
  case("msd")
    action = CALC_MSD
    i=i+1  ; call check_arg(i, "INPUT_FILE")
    call get_command_argument(i,infile)
    i=i+1  ; call check_arg(i, "OUTPUT_FILE")
    call get_command_argument(i,outfile)
  case("prune")
    action = PRUNE
    i=i+1  ; call check_arg(i, "INPUT_FILE")
    call get_command_argument(i,infile)
    i=i+1  ; call check_arg(i, "OUTPUT_FILE")
    call get_command_argument(i,outfile)
  case("--period")
    i=i+1  ; call check_arg(i, "PERIOD")
    call get_command_argument(i,option)
    read(option, *) period_pr
  case default
    write(*,'(A,A)') "Unrecognized option: ",trim(option)
    call write_usage()
    stop 2
  end select
end do

if(verb) then
  call write_version()
endif

!$ call set_threads_openmp(nthreads_req)

! Not sure why this is done.
call random_seed(size=ni)
allocate(ssseed(ni))
call cpu_time(dump)
ssseed(1:ni)=modulo(int(dump*dump)+34964287,2345789)

if(para_lec_only) then
  stop 0
end if

! Dispatch to the requested tool.
select case(action)
case(RUN_SIMU)
  call init_simulation(the_sim, fich_dist,fich_par,fich_in)
  if(the_sim%state%inp%semigrand) then
    call sgc_doublerun(the_sim)
  else
    call run_simulation_equil(the_sim, the_sim%state%inp%maxcycle_stab)
    call run_simulation(the_sim, the_sim%state%inp%maxcycle_calc)
  end if
  call deallocate_simulation(the_sim)
case(CALC_GDER)
  call compute_gder(infile, outfile, fich_dist, fich_par,dr_crdf, bstep, estep)
case(CONVERT)
  ! Last argument on command line is the output file
  call convert_conf(ninfiles-1, infiles(1:ninfiles-1), infiles(ninfiles), fich_dist, fich_par, iformat, oformat, ibox)
case(RESCALE)
  call rescale_conf(infile, outfile, factor, fich_dist, fich_par)
case(EXTRACT)
  call extract_conf(infile, outfile, radius_extract, partselected,fich_dist, fich_par)
case(DOUBLE)
  call confdouble(infile, outfile,fich_dist, fich_par)
case(TABLES)
  call do_mktables(dr_mktables,fich_dist, fich_par)
case(CALC_BONDORDER)
  if(use_bocfg) then
    call analyze_bondorder(infile, fich_dist, fich_par, ibox, bo_multiconf, fich_bocfg)
  else
    call analyze_bondorder(infile, fich_dist, fich_par, ibox, bo_multiconf)
  end if
case(CALC_MSD)
  call compute_msd(infile, outfile,fich_dist, fich_par, ibox, bstep, estep)
case(PRUNE)
  call prune_conf(fich_dist, fich_par, infile, outfile, period_pr, estep)
case default
  call write_usage()
  stop 1
end select

deallocate(ssseed)
contains

subroutine check_arg(iarg, desc)
  integer, intent(in) :: iarg
  character(*), intent(in) :: desc
  if(iarg > argc) then
    write(*,'(A,A)') "Syntax Error: Missing Argument : ", desc
    write(*,*)
    call write_usage()
    stop 1
  endif
end subroutine

  subroutine write_usage()
    write(*,'(A)') "Usage: polcolmc [action] [options]"
    write(*,*)
    write(*,'(A)') "Actions:"
    write(*,'(A)') "  run                                        (default) Run a simulation"
    write(*,'(A)') "  rdf      INPUT_FILE OUTPUT_FILE            Calculate a rdf on a 'inst' trajectory with an interval"
    write(*,'(A)') "  convert  INPUT_FILE OUTPUT_FILE            Convert a trajectory/configuration to a different format."
    write(*,'(A)') "  rescale  INPUT_FILE OUTPUT_FILE FACTOR     Rescale a 'inst' configuration by FACTOR"
    write(*,'(A)') "  extract  INPUT_FILE OUTPUT_FILE            Extract the particles around the one of number PART"
    write(*,'(A)') "  mktable  BIN                               Initialize and output potential tables of interval BIN"
    write(*,'(A)') "  bondorder INPUT_FILE                       Compute bond-order parameters and determine &
&the composition of crystal phases."
    write(*,'(A)') "  msd      INPUT_FILE OUTPUT_FILE            Calculate mean square displacements on a 'inst' trajectory."
    write(*,*)
    write(*,'(A)') "Options:"
    write(*,'(A)') "  -h                Print this message and quits"
    write(*,'(A)') "  --version         Print version and quits"   
    write(*,'(A)') "  -v                Increase verbosity"
    write(*,'(A)') "  -l                Only read parameters and print them"
    write(*,'(A)') "  -i FILE           Get initial configuration from *FILE* (default: in.inst)"
    write(*,'(A)') "  -d FILE           Get distribution from file *FILE* (default: distrib)"
    write(*,'(A)') "  -p FILE           Get parameters from file *FILE* (default: polcolmc.cfg)"
    write(*,'(A)') "  --dr DR           (rdf) Interval length"
    write(*,'(A)') "  --bstep BSTEP     (rdf,msd) First step in the trajectory"
    write(*,'(A)') "  --estep ESTEP     (rdf,msd) Last step in the trajectory"
    write(*,'(A)') "  --if FORMAT       (convert) Input format (default: inst)"
    write(*,'(A)') "  --of FORMAT       (convert) Output format (default: xyz)"
    write(*,'(A)') "  --box BOX         (convert) Number of the box to consider (default: 1)"
    write(*,'(A)') "  --particle PART   (extract) Number of the selected particle"
    write(*,'(A)') "  --radius FORMAT   (extract) Radius around the selected particle"
    write(*,'(A)') "  --bondorder-cfg BOND_ORDER_CONFIG  Separate configuration file for bond-order analysis"
    
  end subroutine write_usage
  
  subroutine write_version()
    write(output_unit,'(A,", version ",A)') trim(PROGNAME),POLCOLMC_VERSION
#ifdef __COMMIT_HASH__
    write(output_unit,'(2A)') "Commit hash: ", __COMMIT_HASH__
#endif
    write(output_unit,*)
  end subroutine
 
end program pdcolmc


module simulation
  !! Module helping to do a MC simulation as implemented here
  !! (built-in moves include translations, swaps, volume change and moves for the Gibbs ensemble.)
  !! (see [[moves]] to add user-defined movements).
  !!
  !!### Example
  !!
  !! This initialize and run a simulation of distribution file ``mydistro``,
  !! of configuration file ``mysimu.cfg``.
  !! A blank `positions` parameter (normally positions of
  !! particles) means here that the initial positions will be generated at random
  !! instead of a file.
  !!
  !!```fortran
  !!use simulation
  !!type(MCSimulation) :: simu
  !!! Initialize a simulation.
  !!call init_simulation(simu, distribution="mydistro",parameters="mysimu.cfg", positions="") 
  !!
  !!! Run a 10-step simulation after 10 steps of equilibration
  !!call run_simulation(simu, maxsteps=10, maxsteps_equil=10)
  !!```

  use state
  use moves
  use correls
  use spemoves
  use sdeq_comp
  use fluct_dens
  use conftemp
  use zdens
  use chempot
  use swapmap
  use statcell
  use ergo
  implicit none

  character(*), parameter, private :: fichqvals="qvals.dat"
  character(*), parameter, private :: suff_fluctu_denses="_denses.dat"
  character(*), parameter, private :: suff_fluctu_rad_denses="_rad_denses.dat"
  character(*), parameter, private :: suff_fluctu_phi_denses="_phi_denses.dat"
  character(*), parameter, private :: suff_fluctu_cov_denses="_cov_denses.dat"
 
  character(*), parameter, private :: suff_bigsout="_bigseff.dat"
  character(*), parameter, private :: suff_out_inst="_traj.inst"
  character(*), parameter, private :: suff_out_nrj="_nrg.dat"
  character(*), parameter, private :: suff_out_ct = "_ct.dat"
  character(*), parameter, private :: suff_end_conf="_out.inst"
  !character(*), parameter, private :: suff_p="_nrg.out"
  character(*), parameter, private :: suff_s="_seff.dat"
  character(*), parameter, private :: suff_sfull="_seff_full"
  character(*), parameter, private :: suff_bmu="_bmu.dat"
  character(*), parameter, private :: suff_swapmap="_swapmap.dat"
  character(*), parameter, private :: suff_swapmapbox="_swapmapbox.dat"
  character(*), parameter, private :: suff_g="_rdf.dat"
  character(*), parameter, private :: suff_pop="_pop.dat"
  character(*), parameter, private :: suff_allg="_rdf_all"
  character(*), parameter, private :: suff_init ="_init.inst"
  character(*), parameter, private :: suff_smax = "_smax.dat"
  character(*), parameter, private :: suff_nrgergo = "_nrgergo.dat"
  
  integer, private :: out_nrj(max_boxes)
  integer, private :: out_nrgergo(max_boxes)

  type MCSimulation
    !! Derived type representing a simulation.
    type(PcmcState), allocatable :: state !! simulation state
    type(MCMoves) :: moves   !! moves
    character(:), allocatable :: prefix
    logical :: prod
    type(NrgRoutines) :: nrg
    integer :: n_g,n_s
    type(CorrelsRec), allocatable :: res_rdf(:) !! rdf results and temps 
    type(StructqRec), allocatable :: res_sq(:)  !! S(q) results
    real(8), allocatable :: Smax(:)
    real(8), allocatable :: epot_tail(:),virp_tail(:)
    type(NrgErgo), allocatable :: nea(:)
  end type MCSimulation

  private :: run_simulation_nofinaloutput, output_seff, output_sfullallg, output_maxs, output_rdf
  private :: output_nrg, output_pop, output_conf, output_nrg_file, calcul_tail, update_dx_trans
  private :: NrgErgo, nrgergo_update
  private :: allocate_simulation_intern

contains

  subroutine set_threads_openmp(nthreads_req)
    !! Sets the number of threads used by the simulation.
    !$ use omp_lib
    integer, intent(in) :: nthreads_req
    !$ integer :: nthreads
    !$ if(nthreads_req>0) then
    !$   call omp_set_dynamic(.false.)
    !$   call omp_set_num_threads(nthreads_req)
    !$ endif
    !$ nthreads=omp_get_max_threads()

    !$ write(*,"(A,I0,A,I0,A)") "Working with ", nthreads, " thread(s) (",omp_get_num_procs()," CPUs available)."
  end subroutine set_threads_openmp

  subroutine allocate_simulation_intern(sim, maxbox)
    type(MCSimulation) :: sim
    integer, intent(in) :: maxbox
    allocate(sim%state)
    allocate(sim%res_rdf(maxbox), sim%res_sq(maxbox))
    allocate(sim%Smax(maxbox), sim%epot_tail(maxbox), sim%virp_tail(maxbox))
    allocate(sim%nea(maxbox))
  end subroutine

  subroutine deallocate_simulation(sim)
    type(MCSimulation) :: sim
    deallocate(sim%state)
    deallocate(sim%res_rdf, sim%res_sq)
    deallocate(sim%Smax, sim%epot_tail, sim%virp_tail)
    deallocate(sim%nea)
    if(allocated(sim%prefix)) deallocate(sim%prefix)
    call deallocate_moves(sim%moves)
  end subroutine

  subroutine init_modes_nrg_cutoff(st, nr)
    use nrg
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    real(8) ::  dump,het,volume,inv_debye_length
    integer :: i,ntotpar,k

    ntotpar = sum(st%dist%pop(:))

    call ran2_put_seed(st%rng, st%inp%seed)
    call ran2_put_seed(st%rng_an, st%inp%seed_analysis)

    st%rpar%dx_trans=min(st%rpar%dx_trans/st%rpar%Lmin, 0.20D0)
    if(st%mode_box==2) then
      if(st%rpar%dx_trans_wall < 0) then ! default value = -1
        st%rpar%dx_trans_wall = st%rpar%dx_trans
      else
        st%rpar%dx_trans_wall = min(st%rpar%dx_trans_wall/st%rpar%Lmin, 0.20D0)
      endif
    endif
    call nrg_init(st, nr)

    ! Set particles props (ech). Needs to to that after st%ndat%eff_charge is set
    do k=1,st%ntotbox
      do i=1,st%boxes(k)%n
        call part_set_props(st%boxes(k)%parts(i), st)
      end do
    end do
    
    het=0.D0
    do i=1,st%dist%nfam
      het=het+st%dist%pop(i)*(st%dist%rad(i))**3
    enddo
    volume = sum(st%boxes(1:st%ntotbox)%volume)
    het=het*(4*PI/3)/volume

    inv_debye_length=1.D0/st%inp%debye_length
  
    call init_cutoff(st, nr)

  end subroutine init_modes_nrg_cutoff

  subroutine init_cutoff(st, nr)
    use nrg
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    integer :: i
    real(8) :: dump

    do i=floor(st%rpar%Lmin/(2*st%inp%dr_calc_cutoff)),1,-1
      st%rc=i*st%inp%dr_calc_cutoff
      if(abs(nr%potential(st%rc,st%dist%nfam,st%dist%nfam,st))>=st%inp%eps*kt) exit
    enddo
    dump = maxval(st%dist%rad(1:st%dist%nfam))
    if(st%rc<2*dump) then
      st%rc = 2*dump+st%inp%dr_calc_cutoff
    end if
  !  
    if(st%inp%rcutoff>0.D0) then
      st%rcut_ap=st%rc
      st%rc=st%inp%rcutoff
    endif
    st%rc = min(st%rc, st%ndat%rcutoff_max)  ! Needed to handle intrinsic limits (Yukawa2...)

    st%rcrc=st%rc*st%rc
  end subroutine

  subroutine init_positions(st,fichin, unit)
    !! Initialize particle positions
    use random2
    use misc, only: random_placement
    use iso_fortran_env
    implicit none
    type(PcmcState) :: st
    character(*), intent(in) :: fichin
    integer, intent(in) :: unit
    real(8) :: vec3(3),boxtens(3,3),boxstar(3,3)
    real(8) :: dump, volume
    logical :: inestla
    integer :: i,j,k, iostat,Ntotpar,ichunk
    type(Pmconf) :: fin

    ! Can be overriden later by initial conf reading
    st%step = st%inp%initial_step
    inquire(file=fichin,exist=inestla)
    if(.not.inestla) then
      write(unit,'(A,A)') trim(fichin)," cannot be found. Going to compute an initial configuration now..."
      call init_box_from_configuration(st,st%inp%boxLx,st%inp%boxLy,st%inp%boxLz,&
        st%inp%boxalpha,st%inp%boxbeta,st%inp%boxgamma,st%inp%phi_inp, unit)

      ! temporary placement in box 1
      i=0
      do j=1,st%dist%nfam
        do while(i<sum(st%dist%pop(1:j)))
          i=i+1
          st%boxes(1)%Parts(i)%Famille=j
          st%boxes(1)%Parts(i)%rayon=st%dist%rad(j)
          st%boxes(1)%Parts(i)%numero=i
        enddo
      enddo

      ! places particles in the boxes by chunk Nptot / nboxes
      Ntotpar = i
      ichunk=Ntotpar/st%ntotbox
      st%boxes(1)%n=ichunk
      do i=2,st%ntotbox-1
        st%boxes(i)%n=ichunk
        do j=1,ichunk
          st%boxes(i)%Parts(j) = st%boxes(1)%Parts((i-1)*ichunk+j)
        enddo
      end do
      if(st%ntotbox>1) then
        st%boxes(st%ntotbox)%n=Ntotpar-ichunk*(st%ntotbox-1)
        do j=1,Ntotpar-ichunk*(st%ntotbox-1)
          st%boxes(st%ntotbox)%Parts(j) = st%boxes(1)%Parts((st%ntotbox-1)*ichunk+j)
        end do
      endif

      do k=1,st%ntotbox
        if(random_placement(st%boxes(k)%n,st%boxes(k)%Parts,st%inp%seed,st%boxes(k)%met,st%mode_box)) then
          write(unit,'(A)') "...Found."
        else
          write(unit,'(A)') "...Not Found. Error, exiting."
          stop 7
        end if
      end do

      ! Needs to be called here inc case an input file disagrees.
      call affect_radius_from_distrib(st)
    else
      fin = open_pmconf(fichin,old=.true.)
      call read_pmconf_frame(st,fin,.true., iostat)
      ! Takes the last frame in the file
      do while(iostat /= iostat_end)
        call read_pmconf_frame(st,fin,.false., iostat)
        if(iostat /= 0 .and. iostat /= iostat_end) then
          write(error_unit,'(3A,I0,A)') "Error: Problem in '", trim(fichin), "'(code ", iostat, ")."
          stop 6
        end if
      end do
      call close_pmconf(fin)

      ! put everything with positive coordinates in case of
      do k=1,st%ntotbox
        do i=1,st%boxes(k)%n
          do j=1,3
            dump = st%boxes(k)%parts(i)%pos(j)
            st%boxes(k)%parts(i)%pos(j) = modulo(modulo(dump,1.d0),1.d0)
            ! Does not seem to suffice, so:
            if(st%boxes(k)%parts(i)%pos(j)>=  1.d0-1.d-12) then
              st%boxes(k)%parts(i)%pos(j) = 0.d0 
            end if
          end do
        end do
      end do
    end if
    
    do k=1,st%ntotbox
      call calc_metric(st%boxes(k)%tens, st%boxes(k)%met)
      st%boxes(k)%num = k
    end do

    vec3(1)=maxval(abs(st%boxes(1:st%ntotbox)%tens(1,1)))
    vec3(2)=maxval(abs(st%boxes(1:st%ntotbox)%tens(2,2)))
    vec3(3)=maxval(abs(st%boxes(1:st%ntotbox)%tens(2,2)))
    st%Lmax = maxval(vec3)
    vec3(1)=minval(abs(st%boxes(1:st%ntotbox)%tens(1,1)))
    vec3(2)=minval(abs(st%boxes(1:st%ntotbox)%tens(2,2)))
    vec3(3)=minval(abs(st%boxes(1:st%ntotbox)%tens(2,2)))
    st%rpar%Lmin = minval(vec3)

    call allocate_multires(st)

    do k=1,st%ntotbox

      do j=1,st%dist%nfam
        st%boxes(k)%fam(j) = 0
      end do
      do i=1,st%boxes(k)%n
        j=st%boxes(k)%parts(i)%famille
        st%boxes(k)%fam(j) = st%boxes(k)%fam(j) + 1
      end do

    end do

    do k=1,st%ntotbox
      boxtens=st%boxes(k)%tens
      volume = inverse_box(boxtens,boxstar)
      if(volume<=1.D-300) then
        write(error_unit,*) "Error: Box",k,": 'volume<=1.E-300'. Too low value."
        do i=1,3
          write(*,'(3G16.8)') real(boxtens(i,:))
        end do
        stop 6
      end if
      st%boxes(k)%volume = volume
    end do

    write(unit,'(A)') "Volume fractions:"
    do k=1, st%ntotbox
      st%res(k)%hetv=0.0D0
      do j=1,st%dist%nfam
        st%res(k)%hetv=st%res(k)%hetv+st%boxes(k)%fam(j)*st%dist%rad(j)**3
      end do
      st%res(k)%hetv = st%res(k)%hetv * 4.d0*pi/3.d0 /st%boxes(k)%volume
      write(unit,'("Box ", I3, "    ", F15.8)') k, st%res(k)%hetv
    end do

  end subroutine init_positions

  subroutine init_simple_simulation(sim,fichdist,fichpar,fichin)
    !! Initialize a simplified (one run, no analysis) simulation.
    use readparam
    implicit none
    type(MCSimulation) :: sim
    character(*), intent(in) :: fichdist,fichpar,fichin
    logical :: inestla
    integer :: stat
    call allocate_simulation_intern(sim, max_boxes)
      ! Read distribution file
    call init_distribution(sim%state%dist,fichdist, stat)
    if(stat<0) then
      write(*,'(A,A,A)') "Error: distribution '",trim(fichdist),"' cannot be found."
      stop 1
    endif
    write(*,'(A)') "Distribution read."
    sim%state%Nbuff = sum(sim%state%dist%pop) + 10

    call init_input(sim%state%inp)
    ! Read configuration file
    inquire(file=fichpar,exist=inestla)
    if(.not. inestla) then
      write(*,'(A,A,A)') "Error: Configuration file '", trim(fichpar),"' cannot be found."
      stop 1
    end if
    call init_param_list(sim%state%config)
    call read_cfg_file(sim%state,fichpar)
    write(*,'(A)') "Configuration read."
    call configuration_finish(sim%state)

    call init_modes(sim%state)

    call init_positions(sim%state,fichin, output_unit)
    
    call init_modes_nrg_cutoff(sim%state, sim%nrg)
    call add_move(sim%moves, "translations", 1.d0, translation1P)
  end subroutine init_simple_simulation


  logical function check_period(rstep, period)
    integer, intent(in) :: rstep
    integer, intent(in) :: period
    if(period > 0) then
      check_period =  (mod(rstep, period) == 0)
    else
      check_period = .false.
    end if
  end function check_period
  

  subroutine init_step_series(sim,fich_out_inst,fich_out_prefix,first_step)
    !! Initialize a series of MC cycles.
    use iso_fortran_env
    use nrg
    use misc, only: isnot_overlapped
    type(MCSimulation) :: sim !! Simulation
    character(*),intent(in) :: fich_out_inst !! Path to the trajectory file
    character(*),intent(in) :: fich_out_prefix  !! Prefix of the path to the energy files 
    integer, intent(in),optional :: first_step !! Initial step
    real(8) :: dumpv(3),bm(3,3)
    integer :: k,id,jd
    logical :: recouv
    character(512) :: bf
    associate(st => sim%state)
    
    st%sdata%cycle_update_trans = .not. sim%prod ! In production this is not strictly correct
    if(present(first_step)) then
      st%sdata%cycle_first_step = first_step
      st%step = first_step
    else
      st%sdata%cycle_first_step = 0
      st%step = 0
    end if
    st%sdata%nnn = sum(st%boxes(1:st%ntotbox)%n) 
  
    if(st%celllist) then
      do k=1,st%ntotbox
        bm = st%boxes(k)%tens
        dumpv =  [st%rc/bm(1,1), st%rc/bm(2,2), st%rc/bm(3,3)]
        call celldec_cree(st%clistdecomp(k),st%boxes(k)%n,st%boxes(k)%parts,dumpv)
      end do
    end if
  
    do k=1,st%ntotbox
      st%boxes(k)%ener = sim%nrg%energy_box(st,st%boxes(k), st%boxes(k)%press)
      if(.not. st%boxes(k)%ener%noverlap) then
          recouv = isnot_overlapped(k,st%boxes(k)%n,st%boxes(k)%parts,st%boxes(k)%Met,st%rcrc,id, jd,st%mode_box)
          if(.not. recouv) then
  
             write(error_unit,'(A,I0,A,I0,A)') "Fatal Error: step_series_next*() => Overlap on initial position. ("&
                  , id , ",", jd, ")"
             write(error_unit,'("pos(",I0,")=", 3G16.8, " red=", 3G16.8 )') id, &
                  matvec3_prod(st%boxes(k)%tens, st%boxes(k)%parts(id)%pos), st%boxes(k)%parts(id)%pos
             write(error_unit,'("pos(",I0,")=", 3G16.8," red=", 3G16.8 )') jd, &
                  matvec3_prod(st%boxes(k)%tens, st%boxes(k)%parts(jd)%pos), st%boxes(k)%parts(jd)%pos
             if(st%celllist) then
                write(error_unit,'("cell(",I0,")=", G16.8 )') id, st%boxes(k)%parts(id)%cellule
                write(error_unit,'("cell(",I0,")=", G16.8 )') jd, st%boxes(k)%parts(jd)%cellule
             end if
             stop 3
          end if
       end if
    end do
  
    if(st%mode_box == 2) then
      st%transss = 0
      st%reuss_trans_wall = 0
      st%reuss_trans_outwall = 0
      st%trans_outwall = 0
      st%trans_wall = 0
  
      st%sdata%reuss_trans_wall_p=0
      st%sdata%reuss_trans_outwall_p=0
      st%sdata%trans_wall_p=0
      st%sdata%trans_outwall_p=0
    end if
    
    do k=1, st%ntotbox
      write(bf, '(A,"_box",I0,A)') fich_out_prefix, k, suff_out_nrj 
      open(newunit=out_nrj(k),file = bf)
    end do
    call append_nrj(st,out_nrj)
    st%finst = open_pmconf(fich_out_inst, bin=st%inp%use_bin_traj)
    call write_pmconf_frame(st,st%finst)
  
    st%res(1:st%ntotbox)%Ep=st%boxes(1:st%ntotbox)%ener%ep
    st%res(1:st%ntotbox)%hsp=st%boxes(1:st%ntotbox)%press%hsp
    st%res(1:st%ntotbox)%virp=st%boxes(1:st%ntotbox)%press%virp
    
    st%res(1:st%ntotbox)%Var_Ep=st%boxes(1:st%ntotbox)%ener%ep**2
    st%res(1:st%ntotbox)%var_hsp=st%boxes(1:st%ntotbox)%press%hsp**2
    st%res(1:st%ntotbox)%var_virp=st%boxes(1:st%ntotbox)%press%virp**2
  
    do k=1,st%ntotbox
       st%res(k)%popm(1:st%dist%nfam) = 1.d0*st%boxes(k)%fam(1:st%dist%nfam)
       st%res(k)%volm = st%boxes(k)%volume
       st%res(k)%vol2m = st%boxes(k)%volume**2
    end do
    st%sdata%step_count = 1
    st%sdata%step_nrj_n = st%inp%step_nrj
    st%sdata%step_nrj_count =  1
    st%sdata%step_uptr = st%inp%per_reuss_trans
    st%sdata%step_w_n = st%inp%step_w  
  
    st%sdata%trans_eff=0
    st%sdata%reuss_trans=0
    st%sdata%reuss_trans_p=0
    st%sdata%trans_eff_p=0
  
  
    end associate
  end subroutine init_step_series
  

  subroutine update_trans_ampl(st)
    type(PcmcState) :: st
    integer :: rstep
    rstep = st%step-st%sdata%cycle_first_step
    if(st%sdata%cycle_update_trans .and. check_period(rstep, st%inp%per_reuss_trans)) then
      if(st%mode_box == 2) then
        call update_dx_trans(st,st%reuss_trans_outwall-st%sdata%reuss_trans_outwall_p,&
          st%trans_outwall-st%sdata%trans_outwall_p,st%rpar%dx_trans)
        st%sdata%reuss_trans_outwall_p=st%reuss_trans_outwall
        st%sdata%trans_outwall_p=st%trans_outwall

        call update_dx_trans(st,st%reuss_trans_wall-st%sdata%reuss_trans_wall_p,st%trans_wall-st%sdata%trans_wall_p,&
          st%rpar%dx_trans_wall)
        st%sdata%reuss_trans_wall_p=st%reuss_trans_wall
        st%sdata%trans_wall_p=st%trans_wall
        if(st%rpar%dx_trans_wall>=0.25) st%rpar%dx_trans_wall=0.25
      else
        call update_dx_trans(st,st%sdata%reuss_trans-st%sdata%reuss_trans_p,st%sdata%trans_eff-st%sdata%trans_eff_p,&
          st%rpar%dx_trans)
        st%sdata%reuss_trans_p=st%sdata%reuss_trans
        st%sdata%trans_eff_p=st%sdata%trans_eff
      end if
      st%sdata%step_uptr=st%sdata%step_uptr+st%inp%per_reuss_trans
      if(st%rpar%dx_trans>=0.25) st%rpar%dx_trans=0.25
    endif
  end subroutine

  subroutine update_base_res(sim)
    use iso_fortran_env
    type(MCSimulation) :: sim
    integer :: rstep,i,sel_mvmt,k
    real(8) :: emean, evar
    type(NrgRec) :: ener_dump, ener_ip
    associate(st => sim%state)
    rstep = st%step-st%sdata%cycle_first_step
    st%sdata%step_count = st%sdata%step_count + 1
    do k=1,st%ntotbox
      st%res(k)%ep=st%res(k)%ep+st%boxes(k)%ener%ep
      st%res(k)%var_ep=st%res(k)%var_ep+st%boxes(k)%ener%ep**2

      st%res(k)%popm(1:st%dist%nfam) = st%res(k)%popm(1:st%dist%nfam) + 1.d0*st%boxes(k)%fam(1:st%dist%nfam)
      st%res(k)%volm = st%res(k)%volm + st%boxes(k)%volume
      st%res(k)%vol2m = st%res(k)%vol2m + st%boxes(k)%volume**2
    end do

    if(rstep==st%sdata%step_nrj_n) then

      do k=1,st%ntotbox
        if(st%inp%nrg_ergo_analysis) then
          ! The obvious problem here is that it becomes dubious when particle vectors are reorganized
          ener_dump%ep = 0 ; ener_dump%noverlap=.true.
          do i=1,st%boxes(k)%n
            ener_ip = sim%nrg%energy_1p(st, st%boxes(k), st%boxes(k)%parts(i))
            ener_dump = ener_dump + ener_ip
            if(.not. ener_dump%noverlap) exit
            sim%nea(k)%cumepp(i) = sim%nea(k)%cumepp(i) + ener_ip%ep ! maybe should use inst average
          end do
          ener_dump%ep = 0.5d0*ener_dump%ep
          call nrgergo_update(sim%nea(k), st%boxes(k))
        else
          ener_dump = sim%nrg%energy_box(st,st%boxes(k),st%boxes(k)%press)
        end if
        if(.not. ener_dump%noverlap) then
          write(error_unit,'(A,I0,A)') "FatErr:: step_series_next*() => Overlap step ",st%step, " !"
          stop 666
        end if
          if(abs(ener_dump%ep - st%boxes(k)%ener%ep)>st%inp%delta_nrj_tol*abs(ener_dump%ep)&
          .and. ener_dump%ep /= 0.d0) then
              write(error_unit,'(A, G15.8,A,G15.8,A,I0,A,I0,A)') &
                  "In step_series_next(), Energy check :  Recalculated = ",&
                  ener_dump%ep, " Expected = ",st%boxes(k)%Ener%ep, &
                  " (box ",k,", step ",st%step," )"
              stop 8
          end if
        st%boxes(k)%Ener = Ener_dump
          
        st%res(k)%hsp=st%res(k)%hsp+st%boxes(k)%press%hsp
        st%res(k)%var_hsp=st%res(k)%var_hsp+st%boxes(k)%press%virp**2

        st%res(k)%virp=st%res(k)%virp+st%boxes(k)%press%virp
        st%res(k)%var_virp=st%res(k)%var_virp+st%boxes(k)%press%virp**2
      end do
      st%sdata%step_nrj_count = st%sdata%step_nrj_count + 1
      
          
      call append_nrj(st,out_nrj)

      if(st%inp%nrg_ergo_analysis) then
        do k=1,st%ntotbox
          call nrgergo_write(sim%nea(k), rstep-st%inp%step_nrj, out_nrgergo(k))
        end do
      end if
      st%sdata%step_nrj_n=st%sdata%step_nrj_n+st%inp%step_nrj
    end if
        
    if(rstep==st%sdata%step_w_n) then
        call write_pmconf_frame(st,st%finst)
        st%sdata%step_w_n=st%sdata%step_w_n+st%inp%step_w
        flush(st%finst%unit)
    endif
    end associate

  end subroutine

  subroutine step_series_next(sim,step)
    !! Performs a step (equivalent of *n* translations, where *n* is the total number of
    !! particles).
    !! 
    !! @note [[init_step_series]] should be called once before calls to this function. @endnote
    use moves
    use iso_fortran_env
    type(MCSimulation) :: sim !! Simulation
    integer, intent(in), optional :: step !! step number
    integer :: i, sel_mvmt
    if(present(step)) then 
      sim%state%step = step
    else
      sim%state%step = sim%state%step + 1
    end if
  
    call update_trans_ampl(sim%state)

    do i=1,sim%state%sdata%nnn+1
      sel_mvmt = select_move(sim%moves,ran2(sim%state%rng))
      call new_movinst(sim%state, sim%nrg, sim%moves%mv(sel_mvmt))
    enddo

    call update_base_res(sim)
  end subroutine step_series_next

  subroutine update_dx_trans(st,idump,iidump, dx)
    !! Updates translation amplitudes (for equilibration, not strictly correct otherwise)
    implicit none
    type(PcmcState) :: st
    integer(8),intent(in) :: idump  !! number of accepted translations
    integer(8),intent(in) :: iidump !! total number of translations
    real(8), intent(inout) :: dx    !! amplitude to update
    if(iidump > 0) then
      if((dble(idump))/(dble(iidump)) < st%inp%reuss_trans_ref) then
        Dx=reduit_dx*Dx
      else
        Dx=augment_dx*Dx
      endif
    end if
  end subroutine update_dx_trans

  subroutine step_series_end(sim)
    !! Should be called after a series of calls to \ref step_series_next. 
    type(MCSimulation) :: sim
    integer :: k,j
    associate(st => sim%state)

    do k=1,st%ntotbox
      close(out_nrj(k))
    end do
    call close_pmconf(sim%state%finst)


    st%res(1:st%ntotbox)%Ep= st%res(1:st%ntotbox)%Ep/st%sdata%step_count
    
    st%res(1:st%ntotbox)%Var_Ep=st%res(1:st%ntotbox)%Var_Ep/st%sdata%step_count-&
      st%res(1:st%ntotbox)%Ep**2

    st%res(1:st%ntotbox)%volm = st%res(1:st%ntotbox)%volm/st%sdata%step_count
    st%res(1:st%ntotbox)%vol2m = st%res(1:st%ntotbox)%vol2m/st%sdata%step_count - &
      st%res(1:st%ntotbox)%volm**2
  
    
    st%res(1:st%ntotbox)%hsp=st%res(1:st%ntotbox)%hsp/st%sdata%step_nrj_count
    st%res(1:st%ntotbox)%var_hsp=st%res(1:st%ntotbox)%var_hsp/st%sdata%step_nrj_count-&
      st%res(1:st%ntotbox)%hsp**2

    st%res(1:st%ntotbox)%virp=st%res(1:st%ntotbox)%virp/st%sdata%step_nrj_count
    st%res(1:st%ntotbox)%var_virp=st%res(1:st%ntotbox)%var_virp/st%sdata%step_nrj_count-&
      st%res(1:st%ntotbox)%virp**2

    do k=1,st%ntotbox
      st%res(k)%popm(1:st%dist%nfam) = st%res(k)%popm(1:st%dist%nfam)/st%sdata%step_count
      st%res(k)%nm = sum(st%res(k)%popm(1:st%dist%nfam))
      
      st%res(k)%p = (kt*st%res(k)%nm)/(st%res(k)%volm)*UnitPsT*st%inp%temp+st%res(k)%hsp+st%res(k)%virp
      st%res(k)%e = 1.5D0*kt*st%res(k)%nm + st%res(k)%Ep
      st%res(k)%e2c = 1.5D0*kt*kt*st%res(k)%nm+st%res(k)%Var_Ep
    end do
    do k=1,st%ntotbox
      st%res(k)%hetv=0.0D0
      do j=1,st%dist%nfam
          st%res(k)%hetv=st%res(k)%hetv+st%res(k)%popm(j)*st%dist%rad(j)**3
      end do
      st%res(k)%hetv = st%res(k)%hetv * 4.d0*pi/3.d0 /st%res(k)%volm
    end do
    end associate
  end subroutine step_series_end


  subroutine write_move_statistics(sim, unit)
    !! Prints statictics about mouvements (translation and user-defined).
    use moves
    type(MCSimulation), intent(in) :: sim
    integer, intent(in) :: unit !! unit of the file to append
    integer :: i

    write(unit,'(A,I10,A)') 'Statistics over',sim%state%sdata%step_count-1,' cycles :'
    do i=1,sim%moves%nmoves
      if(sim%moves%mv(i)%comp>0) write(unit,'(A,F6.2,A,F6.2,A,F6.2)') sim%moves%mv(i)%name,(100.D0*sim%moves%mv(i)%acc)&
           /(1.D0*sim%moves%mv(i)%comp),'% accepted,  HC rej : ', (100.D0*sim%moves%mv(i)%hcrej)/(1.D0*sim%moves%mv(i)%comp),&
           " CIR rej : ",(100.D0*sim%moves%mv(i)%cirrej)/(1.D0*sim%moves%mv(i)%comp)
    end do
    if(sim%state%mode_box == 2) then
      write(unit,'(A,F6.2, A,F6.2,A,F6.2,A)') ' Wall/(Bulk+Wall)= ', (1.0*sim%state%trans_wall)/(1.D0*sim%state%transss)&
        , ' ,Accepted trans: Wall= ',&
          (100.D0*sim%state%reuss_trans_wall)/(1.D0*sim%state%trans_wall),'% , Bulk=  ', &
            (100.D0*sim%state%reuss_trans_outwall)/(1.D0*sim%state%trans_outwall)  ,'%'
    end if

  end subroutine write_move_statistics

  subroutine correls_init_all(sim)
    type(MCSimulation) :: sim
    integer :: k
    sim%n_g=floor(0.5D0*sim%state%rpar%Lmin/sim%state%inp%dr_rdf)
    sim%n_s=floor((PI*(1.D0/(8*sim%state%inp%dr_rdf)-2.D0/sim%state%rpar%Lmin))/sim%state%inp%dq)
    do k=1,sim%state%ntotbox
      call correls_init(sim%res_rdf(k),sim%state%dist%nfam,sim%state%inp%dr_rdf,0.D0,0.5D0*sim%state%rpar%Lmin)
      call correls_reset(sim%res_rdf(k))
    end do
  end subroutine correls_init_all
  
  subroutine output_conf(st,fich_conf)
    !!  Outputs state to a file (positions...)
    type(PcmcState) :: st
    character(*), intent(in) :: fich_conf
    type(Pmconf) :: fout
    fout = open_pmconf(fich_conf)
    call write_pmconf_frame(st,fout)
    call close_pmconf(fout)
  end subroutine output_conf

  subroutine output_pop(st,fich_pop)
    !! Outputs mean family populations
    type(PcmcState) :: st
    character(*), intent(in) :: fich_pop
    integer :: i,k, ufo
    open(newunit=ufo,file=fich_pop)
    do i=1,st%dist%nfam
      write(ufo,'(G16.8)',advance='no') st%dist%rad(i)
      do k=1,st%ntotbox
        write(ufo,'(G16.8)',advance='no') st%res(k)%popm(i)
      end do
      write(ufo,*)
    enddo

    close(ufo)
  end subroutine output_pop

  subroutine output_nrg(st,fu,k)
    !! Outputs some averaged quantities (energy, pressure ...) to file unit fu and for box k
    type(PcmcState) :: st
    integer, intent(in) :: fu
    integer :: k
    real(8) :: ue
    ue = 1.d0

    write(fu,'(A30," = ",G16.8)') 'Mean Number of particles',st%res(k)%nm
    write(fu,'(A30," = ",G16.8,A)') 'Mean Volume',(st%res(k)%volm)**(1.d0/3.d0),"^3"
    write(fu,'(A30," = ",G16.8)') 'Volume fraction ϕ' , st%res(k)%hetv  
    write(fu,'(A30," = ",G16.8)') 'Ideal Pressure (Pa)',(kt*st%res(k)%nm)/(st%res(k)%volm)*unitPsT*st%inp%temp
    write(fu,'(A30," = ",G16.8,A,G16.8)') 'Hard Sphere Pressure (Pa)' ,st%res(k)%hsp," σ =",&
      sqrt(st%res(k)%var_hsp/(st%sdata%step_count+1))
    write(fu,'(A30," = ",G16.8,A,G16.8)') 'Virial Pressure (Pa)',st%res(k)%virp," σ =",&
      sqrt(st%res(k)%var_virp/(st%sdata%step_count+1))
    write(fu,'(A30," = ",G16.8)') 'Pressure (Pa)',st%res(k)%p
    write(fu,'(A30," = ",G16.8,A,G16.8)') 'Potential Energy p. part. (k_B T)', st%res(k)%ep*ue/st%res(k)%nm,&
      " σ =",sqrt(st%res(k)%var_ep/(st%sdata%step_count+1))/&
      st%res(k)%Nm
    write(fu,'(A30," = ",G16.8)') 'Energy p. part. (k_B T)',st%res(k)%e*ue/st%res(k)%nm

  end subroutine output_nrg

  subroutine output_nrg_file(st,file)
    !! Prints some averaged quantities (energy, pressure ...) to file file
    type(PcmcState) :: st
    character(*), intent(in) :: file
    integer :: k, ufo
    open(newunit=ufo,file=file)
    do k=1,st%ntotbox
      write(ufo,'(8G15.7)',advance='no') st%res(k)%hetv, (kt*st%res(k)%nm)/(st%res(k)%volm)*UnitPsT*st%inp%temp, &
        st%res(k)%hsp,st%res(k)%virp,st%res(k)%p,&
        st%res(k)%Ep,st%res(k)%e,st%res(k)%e2c

    end do
    write(ufo,*)

    close(ufo)
  end subroutine output_nrg_file

  subroutine append_nrj(st,out_nrj)
    !! Appends current energy and pressure to the files.
    implicit none
    type(PcmcState) :: st
    integer,intent(in) :: out_nrj(:)

    integer :: k
    real(8) :: ep,vir,hs,ue
    ue = 1.d0

    do k=1,st%ntotbox
      ep = st%boxes(k)%ener%ep * ue
      hs = st%boxes(k)%press%hsp
      vir = st%boxes(k)%press%virp
      write(out_nrj(k),"(I12,5ES16.7)") st%step,ep,vir,hs,real(st%boxes(k)%n),st%boxes(k)%volume
    end do

  end subroutine append_nrj

  subroutine log_step(sim, step, log_unit)
    !! Prints current energy, number of particles, volume for step step.
    type(MCSimulation) :: sim
    integer, intent(in) :: step !! step number
    integer, intent(in) :: log_unit
    character(5) :: zone
    integer, dimension(8) :: tvals
    integer :: k
    call date_and_time(zone=zone,values=tvals)
    write(log_unit,*)
    write(log_unit,'(A, I10,", ",i0,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," ",A)') "Step ", &
      step, tvals(1:3),tvals(5:7),zone
    do k=1, sim%state%ntotbox
      write(log_unit,'(A, I0, A)') "Box ",k, " :"
      write(log_unit,'(A,G15.7, A)') "Potential Energy              = ", &
        sim%state%boxes(k)%Ener%Ep, " kT"
      write(log_unit,'(A,G15.7, A)') "Potential Energy per particle = ", &
        sim%state%boxes(k)%Ener%Ep/sim%state%boxes(k)%n, " kT"
      write(log_unit,'(A,I10, A)'  ) "Number of particles           = ", sim%state%boxes(k)%n, ""
      write(log_unit,'(A,G15.7, A)') "Volume^(1/3)                  = ", (sim%state%boxes(k)%volume)**(1./3.), " nm"
    end do
    flush log_unit
  end subroutine log_step


  subroutine output_rdf(sim, fichg)
    !! Prints mean RDF (for each box) to a file.
    type(MCSimulation) :: sim
    character(*), intent(in) :: fichg
    integer :: i,k, ufo
    character(512) :: filename

    do k=1,sim%state%ntotbox
      write(filename, '(A,"_box",I0,A)') sim%prefix, k, suff_g
      open(newunit=ufo,file=filename)
      do i=1,sim%n_g
        write(ufo,'(ES16.7,1x,ES16.7)') (i-0.5D0)*sim%state%inp%dr_rdf, sim%res_rdf(k)%gmoy(i)
      end do
      close(ufo)
    enddo

  end subroutine

  subroutine output_seff(sim, fichs)
    !! Prints S(q) (for each box) to a file. 
    type(MCSimulation) :: sim
    character(*), intent(in) :: fichs
    integer :: i,k, ufo
    character(512) :: filename

    do k=1,sim%state%ntotbox
      write(filename, '(A,"_box",I0,A)') sim%prefix, k, suff_s
      open(newunit=ufo,file=filename)
      do i=1,sim%n_s
        write(ufo,'(ES16.7,1x,ES16.7)') sim%res_sq(1)%q(i),sim%res_sq(k)%S(i)
      end do
      close(ufo)
    enddo
  end subroutine output_seff

  subroutine output_sfullallg(sim, fichsfull,fichallg)
    !! Prints all RDFs and full S(q) to files (one for each box).
    type(MCSimulation) :: sim
    character(*), intent(in) :: fichsfull !! prefix path to full S(q) files
    character(*), intent(in) :: fichallg  !! prefix path to RDF files
    character(256) :: fichallg_f
    integer :: k,i,j, ufo
    do k=1,sim%state%ntotbox
      write(fichallg_f,'(A,"_box",I0,"_sfull.dat")') trim(fichsfull), k
      open(newunit=ufo,file=fichallg_f,form='FORMATTED')
      do i=1,sim%res_sq(k)%nsfull
        write(ufo,'(3ES16.7)') sim%res_sq(k)%qfull(i), sim%res_sq(k)%sfull(i), sim%res_sq(k)%formf(i)
      end do
      close(ufo)
    end do

    do k=1,sim%state%ntotbox
      write(fichallg_f,'(A,"_box",I0,"_rdfall.dat")') trim(fichallg), k
      open(newunit=ufo,file=fichallg_f,form='FORMATTED')
      write(ufo,*) "#",sim%res_rdf(k)%nr,sim%res_rdf(k)%ng,sim%res_rdf(k)%dr
      do i=1,sim%res_rdf(k)%ng
        write(ufo,'(ES16.7)',advance='no') (1.D0*i-0.5D0)*sim%res_rdf(k)%dr
        do j=1,(sim%state%dist%nfam*(sim%state%dist%nfam+1))/2
          write(ufo,'(ES16.7)',advance='no') sim%res_rdf(k)%g(i,j)
        end do
        write(ufo,*)
      end do
      close(ufo)
    end do
  end subroutine output_sfullallg

  subroutine output_maxs(sim, file)
    !! Prints Smax to a file.
    type(MCSimulation) :: sim
    character(*), intent(in) :: file
    integer :: k, ufo
    open(newunit=ufo,file=file)
    do k=1,sim%state%ntotbox
      write(ufo,'(2ES16.7)',advance='no') sim%state%res(k)%hetv, maxval(sim%res_sq(k)%S)
    end do
    write(ufo,*)

    close(ufo)
  end subroutine output_maxs

  subroutine output_sim_base(sim)
    !! Prints average energies, pressure....
    use iso_fortran_env
    type(MCSimulation) :: sim
    integer :: i,j,k
    associate(st => sim%state)
    if(sim%prod) then
      write(output_unit,'(A)') "Box        NRJ tail         Vir tail"  
      do k=1,st%ntotbox
        write(output_unit,'(I10,G16.8,G16.8)') k,sim%epot_tail(k),sim%VIRP_tail(k)
      end do
      write(output_unit,*)
    end if
    do k=1,st%ntotbox
      write(output_unit,'("Box ",I0)') k
      call output_nrg(st,output_unit,k)
      if(st%inp%prod .and. st%inp%step_job > 0) write(output_unit,'(A30," = ",G16.8)') 'Maximum of the S(q)',sim%smax(k)
    end do
    write(output_unit,*)
    call output_conf(st, sim%prefix // suff_end_conf)
    !call output_nrg_file(st, sim%prefix // suff_p)
    call output_pop(st, sim%prefix // suff_pop)
    end associate
  end subroutine output_sim_base

  subroutine init_simulation(sim, distribution, parameters, positions)
    !! Initialize the simulation.
    !!
    !! Read the distribution file, the parameter file,
    !! the initial  configuration (or finds one if the file does not exist,
    !! and performs some other initialization.
    !!
    !! Must be called before [[run_simulation]] (and much anything else).
    use readparam
    use iso_fortran_env
    implicit none
    type(MCSimulation) :: sim
    character(*), intent(in), optional :: distribution !! Path to a distribution file
    character(*), intent(in), optional :: parameters   !! Path to a parameter file
    character(*), intent(in), optional :: positions    !! Path to a file with an initial particle configuration
    ! logical, intent(in), optional :: verb
    character(len=:), allocatable :: fich_dist_,fich_par_,fich_in_, fich_log_
    integer :: k, ntotpar, stat
    logical :: fexists
    type(Pmconf) :: finit
    call allocate_simulation_intern(sim, max_boxes)
    associate(st => sim%state)

    if(present(distribution)) then
      fich_dist_ = distribution
    else
      fich_dist_ = "distrib"
    end if

    if(present(parameters)) then
      fich_par_ = parameters
    else
      fich_par_ = "polcolmc.cfg"
    end if

    if(present(positions)) then
      fich_in_ = positions
    else
      fich_in_ = "in.inst"
    end if
    
    write(*,'(A)') "Initializing."

    ! Read distribution file
    call init_distribution(st%dist,fich_dist_, stat)
    if(stat<0) then
      write(*,'(A,A,A)') "Error: distribution '",trim(fich_dist_),"' cannot be found."
      stop 1
    endif
    write(*,'(A)') "Distribution read."
    ntotpar = sum(st%dist%pop)
    st%Nbuff = ntotpar + 10

    call init_input(st%inp)
    inquire(file=fich_par_, exist=fexists)
    if(.not. fexists) then
      write(*,'(A,A,A)') "Error: distribution '",trim(fich_par_),"' cannot be found."
      stop 1
    end if
    call read_cfg_file(st, fich_par_)
    call init_modes(sim%state)

    call init_positions(st,fich_in_, output_unit)
    call init_modes_nrg_cutoff(st, sim%nrg)

    sim%prefix = st%inp%file_prefix

    sim%prod = st%inp%prod

    ! register default moves
    call init_moves(sim%state, sim%moves)

    finit = open_pmconf(sim%prefix // suff_init)
    call write_pmconf_frame(sim%state,finit)
    call close_pmconf(finit)

    do k=1,sim%state%ntotbox
      write(output_unit,'(A,3G16.7)') "Box lengths :",sim%state%boxes(k)%tens(1,1),sim%state%boxes(k)%tens(2,2),&
           sim%state%boxes(k)%tens(3,3)
    end do
    write(output_unit,'(A,I0,A,I0)') "Number of particles: ",sum(sim%state%boxes(1:sim%state%ntotbox)%n)&
         , " Number of Families: ",sim%state%dist%nfam

    if(sim%state%rcut_ap>0.D0) then
      write(output_unit,'(A,G16.8)') "A priori Cutoff :",sim%state%rcut_ap
    endif
    write(output_unit,'(A,G16.8)') "Cutoff :",sim%state%rc

    if(st%inp%rvois_impur_rosbth<=0.D0) then
      st%inp%rvois_impur_rosbth = 1.5D0*(sum(st%boxes(1:st%ntotbox)%volume)/ntotpar)**(1.D0/3.D0)
    end if
    ! Protection
    if(.not.st%inp%std_gibbs_volex .and. st%inp%dv_volex>=dv_volex_max) then
      write(error_unit,*) "Error: Parameter 'exchange deltavolume' is too high for volume changes."
      stop 6
    end if

    if(st%inp%n_rospechmax < 0) then
      st%inp%n_rospechmax = st%dist%nfam + 1
    end if
    st%inp%n_rospechmax = max(st%inp%n_rospechmax,1)

    
    write(output_unit, '(A)') 'Initialisation done'
    write(output_unit, *)
    flush(output_unit)
    end associate
  end subroutine init_simulation

  subroutine simulation_analysis_init(sim)
    !! Initialize analysis modules
    type(MCSimulation) :: sim
    integer :: k
    character(512) :: bf
    associate(st => sim%state)
    if(st%inp%statcell) call statcell_init(st)
    if(st%inp%prob_swap > 0.d0) call swapmap_init(st)
    if(st%inp%prob_swap_box > 0.d0) call swapmap_box_init(st)
    if(st%inp%prob_swap_impur > 0.d0) call swapmap_impur_init(st)
    if(st%inp%nrg_ergo_analysis) then
      call nrgergo_init(sim%nea, st)
      do k=1,st%ntotbox
        write(bf, '(A,"_box",I0,A)') sim%prefix,k,suff_nrgergo
        open(newunit=out_nrgergo(k), file=bf)
      end do
    end if

    if(sim%prod) then
      if(st%inp%step_job > 0) then
        call correls_init_all(sim)
        do k=1,st%ntotbox
          call correls_new_step(sim%res_rdf(k),st%boxes(k)%met,st%boxes(k)%n,st%boxes(k)%Parts,st%mode_box)
        end do
      end if
      if(st%inp%step_bigs>0) then
        call sdeq_comp_reset()
        call sdeq_comp_new_step(st%boxes(1)%Parts,st%boxes(1)%n)
      end if

      if(st%inp%chempwid_doit .or. st%inp%widom_period > 0) call chempwid_init(st)

      if(st%inp%pZdens > 0) call zdens_init(st)

      if(st%inp%fluctudens) then
        call fluct_dens_init(st, st%inp%fich_fluctu)
        call fluct_dens_reset()
      end if
      if(st%inp%conf_temp_on) call conftemp_init(st, sim%nrg, sim%prefix // suff_out_ct, st%sdata%cycle_first_step)
    else
      if(st%inp%conf_temp_on) call conftemp_init(st, sim%nrg, sim%prefix // suff_out_ct, st%sdata%cycle_first_step)

    end if
    end associate
  end subroutine simulation_analysis_init

  subroutine simulation_analysis_step(sim)
    !! Makes new step of analysis
    type(MCSimulation) :: sim
    integer :: k, rstep
    associate(st => sim%state)
    rstep = sim%state%step - sim%state%sdata%cycle_first_step
    if(st%inp%statcell .and. check_period(rstep, st%inp%step_nrj)) &
         call statcell_update(st)

    if(sim%prod) then
      if(check_period(rstep, st%inp%step_job)) then
        do k=1,st%ntotbox
          call correls_new_step(sim%res_rdf(k),st%boxes(k)%met,st%boxes(k)%n,st%boxes(k)%Parts,st%mode_box)
        end do
      end if
      if(check_period(rstep, st%inp%step_bigs)) &
           call sdeq_comp_new_step(st%boxes(1)%Parts,st%boxes(1)%n)
      if(check_period(rstep, st%inp%pZdens)) &
           call zdens_update(st)
      if(check_period(rstep, st%inp%stepfluctu)) &
           call fluct_dens_update(st)
      if(st%inp%widom_period > 0) then
        call chempwid_update(st, sim%nrg, st%inp%widom_samples)
      end if
    end if
    if(st%inp%conf_temp_on .and. check_period(rstep, st%inp%step_nrj)) &
         call conftemp_new_step(st, sim%nrg, st%step)

    end associate
  end subroutine simulation_analysis_step

  subroutine simulation_analysis_finish(sim)
    !! Does the finish step of analysis
    use iso_fortran_env
    type(MCSimulation) :: sim
    type(PcmcState) :: st
    integer :: k, itemp
    real(8) :: force2(max_boxes), vseconde(max_boxes) ! not really used atm
    associate(st => sim%state)

    if(st%inp%statcell) then
      call statcell_finish(st)
      call statcell_log(st, output_unit)
    end if
    if(st%inp%prob_swap > 0.d0) call swapmap_finish(st)
    if(st%inp%prob_swap_box > 0.d0) call swapmap_box_finish(st)
    if(n_impur_stat_swap>0) then
      call swapmap_impur_log(st, output_unit)
      write(output_unit, *)
    end if
    if(st%inp%nrg_ergo_analysis) then
      do k=1,st%ntotbox
        close(out_nrgergo(k))
      end do
    end if

    if(sim%prod) then
      if(st%inp%step_job > 0) then
        do k=1,st%ntotbox
          call correls_end_step(sim%res_rdf(k),st%res(k)%volm,st%res(k)%Popm)
          itemp=floor(st%rc/st%inp%dr_rdf)+1
          call calcul_tail(st,sim%nrg,sim%n_g,sim%res_rdf(k)%g,st%dist%nfam,st%inp%dr_rdf,st%res(k)%volm,st%res(k)%popm,itemp,&
               3*itemp,sim%virp_tail(k),sim%epot_tail(k))
          sim%epot_tail(k) = sim%epot_tail(k)
          sim%virp_tail(k) = sim%virp_tail(k)*UnitPsT*st%inp%temp
        end do
        do k=1,st%ntotbox
          if( sim%res_rdf(k)%step_calc > 0) then
            call correls_reduit_g(sim%res_rdf(k),st%res(k)%popm)
          end if
        end do
        do k=1,st%ntotbox
          st%res(k)%virp=st%res(k)%virp+sim%virp_tail(k)
          st%res(k)%P=st%res(k)%P + sim%virp_tail(k)
          st%res(k)%Ep=st%res(k)%Ep + sim%epot_tail(k)
          st%res(k)%E=st%res(k)%E + sim%epot_tail(k)
        end do
        ! exp-like  effective structure factor
        do k=1,st%ntotbox
          call SofqfromGder(sim%res_sq(k),sim%res_rdf(k),st%res(k)%volm,st%res(k)%Popm,&
               st%dist%rad,2*PI/st%rpar%Lmin,PI/(8*st%inp%dr_rdf),st%inp%dq)
          sim%smax(k) = maxval(sim%res_sq(k)%S)
        end do
      end if
      if(st%inp%step_bigs>0) call sdeq_comp_end_step(st%boxes(1)%Parts,st%boxes(1)%n)
      if(st%inp%pZdens > 0) call zdens_finish(st)
      if(st%inp%chempwid_doit .or. st%inp%widom_period > 0) call chempwid_finish(st)
      if(st%inp%fluctudens) call fluct_dens_finish()
    else
      !
    end if
    if(st%inp%conf_temp_on) call conftemp_finish(st%ntotbox,force2, vseconde)
    end associate
  end subroutine simulation_analysis_finish

  subroutine simulation_analysis_output(sim)
    !! Output results from analysis modules
    use iso_fortran_env
    type(MCSimulation) :: sim
    integer :: k
    type(PcmcState) :: st
    associate(st => sim%state)

    if(sim%prod) then
      if(st%inp%step_job > 0) then
        call output_rdf(sim, sim%prefix // suff_g)
        call output_seff(sim, sim%prefix // suff_s)
        call output_sfullallg(sim, sim%prefix, sim%prefix)
        !call output_maxs(sim, sim%prefix // suff_smax)

        if(st%inp%fluctudens) then
          call fluct_dens_output(sim%prefix // suff_fluctu_denses, sim%prefix // suff_fluctu_rad_denses, &
               sim%prefix // suff_fluctu_phi_denses, sim%prefix // suff_fluctu_cov_denses)
        end if
        
        if(st%inp%fluctudens) then
          call fluct_dens_log_results(st, output_unit)
        end if
        
        if(st%inp%pZdens > 0) then
          call zdens_output(st, sim%prefix)
        end if
      end if
  
      if(st%inp%step_bigs>0) call sdeq_comp_output(sim%prefix // suff_bigsout)

      if(st%inp%chempwid_doit .or. st%inp%widom_period > 0) then
        call chempwid_output(st, sim%prefix // suff_bmu)
      end if
    end if
    if(st%inp%prob_swap > 0.d0) then
      call swapmap_output(st, sim%prefix // suff_swapmap)
    end if

    if(st%inp%prob_swap_box > 0.d0) then
      call swapmap_box_output(st, sim%prefix // suff_swapmapbox)
    end if
    end associate
  end subroutine simulation_analysis_output
  
  subroutine simulation_analysis_free(sim)
    !! Free ressources of analysis modules
    type(MCSimulation) :: sim
    integer :: k
    type(PcmcState) :: st
    associate(st => sim%state)
    if(st%inp%statcell) call statcell_free()
    if(st%inp%prob_swap > 0.d0) call swapmap_free()
    if(st%inp%prob_swap_box > 0.d0) call swapmap_box_free()
    if(n_impur_stat_swap>0) call swapmap_impur_free()
    if(st%inp%nrg_ergo_analysis) then
      do k=1,st%ntotbox
        call nrgergo_free(sim%nea(k))
      end do
    end if
    if(sim%prod) then
      if(st%inp%step_job > 0) then
        do k=1,st%ntotbox
          call correls_free(sim%res_rdf(k))
        end do

        do k=1,st%ntotbox
          call structq_free(sim%res_sq(k))
        end do
      end if
      if(st%inp%fluctudens) call fluct_dens_free()
      if(st%inp%pZdens > 0) call zdens_free()
      if(st%inp%chempwid_doit) call chempwid_free()
    end if
    end associate
  end subroutine simulation_analysis_free

  
  subroutine simulation_main_output(sim)
    !! Prints main results from simulation to a file.
    type(MCSimulation) :: sim
    integer :: ofu
    integer :: k,i
    real(8) :: nm, vm,T, ue
    
    open(newunit=ofu, file=sim%prefix // "_res.out")
    write(ofu, '("steps = ",I0)') sim%state%sdata%step_count-1

    ! Move stats
    do i=1,sim%moves%nmoves
      write(ofu, '("[[moves]]")')
      write(ofu, '("name = ",A1,A,A1)') '"',trim(sim%moves%mv(i)%name),'"'
      write(ofu, '("done = ", I0)') sim%moves%mv(i)%comp
      write(ofu, '("accepted = ", ES16.7)') divideOrZero(sim%moves%mv(i)%acc,sim%moves%mv(i)%comp)
      write(ofu, '("hard_core_rej = ", ES16.7)') divideOrZero(sim%moves%mv(i)%hcrej,sim%moves%mv(i)%comp)
      write(ofu, '("cir_rej = ", ES16.7)') divideOrZero(sim%moves%mv(i)%cirrej,sim%moves%mv(i)%comp)
    end do
    write(ofu,*)

    ! Energies etc...
    do k=1, sim%state%ntotbox
      write(ofu, '("[[results]]")')
      write(ofu, '("box = ", I0)') k

      ue = UnitEsT*sim%state%inp%temp
      nm = sim%state%res(k)%nm
      vm = sim%state%res(k)%volm
      T = sim%state%inp%temp
      write(ofu, '("mean_number = ", ES16.7)') nm
      write(ofu, '("mean_volume = ", ES16.7)') vm
      write(ofu, '("volume_fraction = ", ES16.7)') sim%state%res(k)%hetv
      write(ofu, '("ideal_pressure = ", ES16.7)') nm/vm*unitPsT*T
      write(ofu, '("hard_sphere_pressure = ", ES16.7)') sim%state%res(k)%hsp
      write(ofu, '("hard_sphere_pressure_variance = ", ES16.7)') sim%state%res(k)%var_hsp
      write(ofu, '("virial_pressure = ", ES16.7)') sim%state%res(k)%virp
      write(ofu, '("virial_pressure_variance = ", ES16.7)') sim%state%res(k)%var_virp
      write(ofu, '("pressure = ", ES16.7)') sim%state%res(k)%p
      write(ofu, '("potential_energy = ", ES16.7)') sim%state%res(k)%ep*ue
      write(ofu, '("energy = ", ES16.7)') sim%state%res(k)%e*ue
      write(ofu, '("energy_per_particle = ", ES16.7)') sim%state%res(k)%e*ue/nm
      write(ofu, '("energy_variance = ", ES16.7)') sim%state%res(k)%e2c
      if(sim%prod) then
        write(ofu, '("tail_energy = ", ES16.7)') sim%epot_tail(k)
        write(ofu, '("tail_pressure = ", ES16.7)') sim%virp_tail(k)
        write(ofu, '("Smax = ", ES16.7)') sim%smax(k)
      end if
    end do
    close(ofu)
  contains
    real(8) function divideOrZero(a,b)
      integer(8), intent(in) :: a,b
      if(b /= 0) then
        divideOrZero = (1.0*a)/(1.0*b)
      else
        divideOrZero = 0.0
      end if
      end function

  end subroutine simulation_main_output

  subroutine run_simulation_nofinaloutput(sim, maxsteps)
    !! Runs a simulation (helper for [[run_simulation]]
    use iso_fortran_env
    type(MCSimulation) :: sim
    integer, intent(in) :: maxsteps
    type(PcmcState) :: st
    integer :: first_step, step, rstep, log_unit
    associate(st => sim%state)
    open(newunit=log_unit, file=sim%prefix // "_log.txt")
    call reset_all_moves(sim%moves) ! Maybe not necessary
    first_step = sim%state%step
    call init_step_series(sim,sim%prefix // suff_out_inst, sim%prefix,&
         first_step=first_step)
    call simulation_analysis_init(sim)
    call log_step(sim,first_step, log_unit)
    do step = 1 + first_step, first_step + maxsteps
      call step_series_next(sim,step)
      rstep = step-first_step
      if(check_period(rstep, st%inp%step_log)) &
           call log_step(sim, step, log_unit)
      call simulation_analysis_step(sim)
    end do
    call step_series_end(sim)
    call write_move_statistics(sim, output_unit)
    call simulation_analysis_finish(sim)
    call simulation_main_output(sim)
    call simulation_analysis_output(sim)
    call simulation_analysis_free(sim)
    close(log_unit)
    end associate
  end subroutine run_simulation_nofinaloutput
  
  subroutine run_simulation_equil(sim, maxsteps)
    !! Run a pre_equilibration run of `steps` steps.
    use iso_fortran_env
    type(MCSimulation) :: sim
    integer, intent(in) :: maxsteps
    character(:), allocatable :: prefix
    logical :: prod
    prod = sim%prod
    prefix = sim%prefix
    if(maxsteps > 0) then
      sim%prod = .false.
      sim%prefix = prefix // "_eq"
      call run_simulation_nofinaloutput(sim, maxsteps) 
      write(output_unit, *)
      write(output_unit, '(A)') "Pre-equilibration done"
      write(output_unit, *)
      ! restore val for main run
      sim%prod = prod
      sim%prefix =prefix
      end if
  end subroutine

  subroutine run_simulation(sim, maxsteps)
    !! Runs a simulation run of `maxsteps` steps.
    use iso_fortran_env
    implicit none
    type(MCSimulation) :: sim
    integer, intent(in) :: maxsteps

    if(maxsteps > 0) then
      call run_simulation_nofinaloutput(sim, maxsteps)
      write(output_unit, *)
      write(output_unit, '(A)') "Main run done"

      call output_sim_base(sim)
      write(output_unit, '(A)') "Output Done"
    endif

  end subroutine run_simulation

  subroutine calcul_tail(st,nr,n_g,gr,n_r,dr,volume,pop,nmin,nmax,virp_tail,epot_tail)
    !! Calculates the tail contribution to the energy (`virp_tail`) and pressure (`epot_tail`),
    !! from the RDFs.
    implicit none
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    integer,intent(in) :: n_g,n_r,nmin,nmax
    real(8),intent(in) :: gr,dr,volume,pop
    real(8),intent(inout) :: virp_tail,epot_tail
    dimension :: pop(*),gr(n_g,*)
    integer :: i,j,k,ind
    real(8) :: r, fac
    epot_tail=0.d0
    virp_tail=0.d0
    do k=nmin,nmax
      r=(k-0.5d0)*dr
      do i=1,n_r ; do j=i,n_r
        if(i==j) then
          ind=i
          fac = 0.5d0
        else
          ind=((2*n_r-i-1)*i)/2+j
          fac = 1.0
        endif
        if(k<=n_g) then
          epot_tail=epot_tail+fac*4*pi*dr*r*r*gr(k,ind)*pop(i)*pop(j)*nr%potential(r,i,j,st)/volume
          virp_tail=virp_tail+fac*4*pi*dr*r*r*r*gr(k,ind)*pop(i)*pop(j)*nr%force(r,i,j,st)/(3*(volume**2))
        else
          epot_tail=epot_tail+fac*4*pi*dr*r*r*pop(i)*pop(j)*nr%potential(r,i,j,st)/volume
          virp_tail=virp_tail+fac*4*pi*dr*r*r*r*pop(i)*pop(j)*nr%force(r,i,j,st)/(3*(volume**2))
        endif
      enddo; enddo
    enddo
  end subroutine calcul_tail

end module simulation

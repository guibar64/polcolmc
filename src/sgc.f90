module sgc
  !! SemiGrand Canonical ensemble: N, T, V or P fixed + Δμ(R) fixed
  
  use ieee_arithmetic
  use simulation
  use readparam
  use histograms
  implicit none

  !! Input parameters for SGC
  type SgcParams
    real(8) :: particle_resize_cvf_prob = 0.0
    real(8) :: particle_resize_prob = 0.0_8
    logical :: big_mvmt = .true.
    integer :: kappa_update_period = 100
    logical :: update_kappa_radius_resize = .false.
    logical :: npt_special = .false.
    real(8) ::  min_radius = 2., max_radius = 30., dr_resize_max = 1.0
    real(8) ::  min_radius_psd = 2., max_radius_psd = 30., dr_psd = 0.10_8
    logical :: update_chempot_table = .false.
    real(8) :: chempot_update_factor = 0.5
  end type

  real(8), allocatable :: kappa_table(:)  ! effective kappas
  real(8), allocatable :: radius_table(:) ! radii for R_cell interpolation
  real(8) :: radius_table_min, radius_table_max
  real(8), allocatable :: rcell3_table(:,:) ! R_cell(R, kappa)
  real(8), allocatable :: rcell3_der2_table(:,:)
  integer, parameter :: sigma_kappa_n = 4, sigma_rad_n = 3
  real(8) :: sigma_coeffs(sigma_kappa_n,sigma_rad_n) ! coeffs of polynomial fit of eff charge

  type(HistoR8), private :: psd, psd_interm
  logical, private :: do_psd=.false.

  type(SgcParams), private :: sgcpars

  type ChemPotData
    real(8), allocatable :: radius(:), table(:), table_der2(:)
  end type

  type(ChemPotData), private :: chemdat

contains

subroutine sgc_doublerun(simu)
  use iso_fortran_env
  type(MCsimulation) :: simu
  logical :: oldprod
  character(len=:), allocatable :: prefix
  real(8) :: prob_volmod

  if(simu%state%ntotbox > 1) then
    write(output_unit, '(A)') "This simulation is not compatible with multiple boxes. Exiting"
    stop 1
  end if

  call treat_additional_config(simu, sgcpars)

  call add_move(simu%moves, "Particle resize" , sgcpars%particle_resize_prob, particle_resize)
  call add_move(simu%moves, "Particle resize (cte Vₚ)" , sgcpars%particle_resize_cvf_prob, particle_resize_cvf)

  if(sgcpars%npt_special) then
    if(simu%state%inp%prob_volmod<=0._8) then
      prob_volmod=1._8/dble(maxval(simu%state%boxes(1:simu%state%ntotbox)%n))
    else
      prob_volmod = simu%state%inp%prob_volmod
    endif
    call add_move(simu%moves, "Volume change (spe)" , prob_volmod ,changevol_npt_sgc)
  end if

  call init_tables("table_fits_keffs.dat")

  ! Initialize charges and kappa
  call update_kappa_from_dist(simu%state, 1)

  ! charges may have changed, need potentially recompute cutoff radius.
  block
    integer :: i, j
    associate(st => simu%state)
    do i=1,st%boxes(1)%n
      j = st%boxes(1)%parts(i)%famille
      st%ndat%eff_charge(j) = st%boxes(1)%parts(i)%ech
    end do
    call init_cutoff(st, simu%nrg)
    write(output_unit, '("Really, cutoff radius = ", G16.7)') st%rc
    end associate
  end block

  block
    real(8) :: kappa, r
    integer :: j, i, fu, nfam
    real(8), allocatable :: rads(:), charges(:), pop(:)
    associate(st => simu%state)
    nfam = st%dist%nfam
    allocate(rads(nfam), charges(nfam), pop(nfam))
    kappa = get_kappa(st,1)
    pop = 0
    do i=1,st%boxes(1)%n
      j = st%boxes(1)%parts(i)%famille
      r = st%boxes(1)%parts(i)%rayon
      rads(j) = r
      charges(j) = st%boxes(1)%parts(i)%ech * exp(- kappa * r)*(1 + kappa*r)
      pop(j) = pop(j) + 1
    end do
    open(newunit=fu, file ="distrib-fromtable")
    write(fu, *) nfam
    do i=1,nfam
      write(fu, '(I10, 2ES16.7)') int(pop(i)), rads(i), charges(i)
    end do
    end associate
  end block


  block
    integer :: npsd, i
    npsd = floor((sgcpars%max_radius_psd - sgcpars%min_radius_psd) / sgcpars%dr_psd) - 1
    call histor8_init(psd, sgcpars%dr_psd, sgcpars%min_radius_psd, npsd)
    do i=1, simu%state%boxes(1)%n
      call histor8_update(psd, simu%state%boxes(1)%parts(i)%rayon)
    end do
 

    call histor8_init(psd_interm, sgcpars%dr_psd, sgcpars%min_radius_psd, npsd)
    do i=1, simu%state%boxes(1)%n
      call histor8_update(psd_interm, simu%state%boxes(1)%parts(i)%rayon)
    end do

  end block

  call init_chempot_tables(chemdat, "chempot_table.dat")

  ! Comp. PSD for prod run only
  do_psd=.false.

  oldprod = simu%prod
  simu%prod = .false.
  prefix = simu%prefix
  simu%prefix = prefix // "_eq"
  call sgc_run(simu, simu%state%inp%maxcycle_stab)

  do_psd = .true.
  simu%prod = oldprod
  simu%prefix = prefix
  call sgc_run(simu, simu%state%inp%maxcycle_calc)
  call output_sim_base(simu)
  call output_psd(psd, simu%prefix // "_meandist.dat", simu%state%res(1)%nm)

  block
    integer :: i
    call histor8_reset(psd)
    do i=1, simu%state%boxes(1)%n
      call histor8_update(psd, simu%state%boxes(1)%parts(i)%rayon)
    end do

    call histor8_reset(psd_interm)
    do i=1, simu%state%boxes(1)%n
      call histor8_update(psd_interm, simu%state%boxes(1)%parts(i)%rayon)
    end do

  end block
  if(sgcpars%update_chempot_table) call update_chempot(psd, chemdat, "distrib_goal.dat", "chempot_table.dat")
end subroutine


subroutine sgc_log_step(sim, step, log_unit)
  type(MCSimulation) :: sim
  integer, intent(in) :: step, log_unit
  call log_step(sim, step, log_unit)
  associate(st => sim%state)

  write(log_unit,'(A,G15.7, A)') "kappa                        = ", get_kappa(st,1), " nm⁻¹"

  end associate
end subroutine

subroutine sgc_run(sim, maxsteps)
  !! Runs a simulation
  type(MCSimulation) :: sim
  integer, intent(in) :: maxsteps
  integer :: first_step, step, rstep, log_unit
  integer :: sel_mvmt, ins, inner_count, i, big_count, big_accepted, smallsteps_count
  integer :: ininnerc, ukap, upsdi
  type(Box) :: oldbox

  !open(newunit=log_unit, file=sim%prefix // "_log.txt")
  log_unit = output_unit
  associate(st => sim%state, nr => sim%nrg)
  call reset_all_moves(sim%moves) ! Maybe not necessary
  first_step = sim%state%step
  call init_step_series(sim,sim%prefix // "_traj.inst", sim%prefix,&
       first_step=first_step)
  call simulation_analysis_init(sim)
  call sgc_log_step(sim, first_step, output_unit)
 
  open(newunit=ukap, file = sim%prefix // "_kappa.dat")
  open(newunit=upsdi, file = sim%prefix // "_intermpsd.dat")
  if(st%inp%internal_steps > 0) then
    st%sdata%nnn = st%inp%internal_steps
  else
    st%sdata%nnn = st%boxes(1)%n
  end if

  if(sgcpars%kappa_update_period > 0) then
    inner_count = sim%state%sdata%nnn / sgcpars%kappa_update_period
    smallsteps_count = sgcpars%kappa_update_period
  else
    inner_count = sim%state%sdata%nnn
    smallsteps_count = 1
  end if
  ! print *, "internal_steps = ", sim%state%sdata%nnn
  big_count = 0
  big_accepted = 0
  oldbox = st%boxes(1) !
  do step = 1 + first_step, first_step + maxsteps
    sim%state%step = step ! 🤔
    call update_trans_ampl(sim%state)
    do ins=1, inner_count
      if(sgcpars%kappa_update_period > 0 .and. ins == inner_count) then
        ininnerc = modulo(sim%state%sdata%nnn, sgcpars%kappa_update_period)
      else
        ininnerc = smallsteps_count
      end if
      do i = 1, smallsteps_count
        ininnerc = ininnerc + 1
        sel_mvmt = select_move(sim%moves, ran2(sim%state%rng))
        call new_movinst(sim%state, sim%nrg, sim%moves%mv(sel_mvmt))
      end do

      if(sgcpars%kappa_update_period > 0) then
        if(sgcpars%big_mvmt) then
          ! κ should change for every radius change, all the interactions in the box would then need
          ! to be recomputed. Instead, the actual Markov chain is done in big steps.
          ! To recover detailed balance, the acceptance ratio between the new and old configuration
          ! is P'eq(o)/P'eq(n)*Peq(n)/Peq(o), which reduce here to Peq(n)/P'eq(n) = exp -β(En - En')
          block
            type(NRGRec) :: epn, en
            real(8) :: oldkp
            integer :: ib = 1

            big_count = big_count + 1
            oldkp = get_kappa(st,1)
            epn =  st%boxes(ib)%ener
            call update_kappa_from_dist(st, ib) ! ⚠ cutoff changes ?
            en = nr%energy_box(st, st%boxes(ib))
            if(en%noverlap .and. metropolis((en%ep-epn%ep)/kT, st%rng)) then
              big_accepted = big_accepted + 1
              st%boxes(ib)%ener = en
              oldbox = st%boxes(ib)
            else
              call set_kappa(st, ib, oldkp)
              st%boxes(ib) = oldbox
              if(st%celllist) call celldec_update_all(st%clistdecomp(ib), st%boxes(ib)%n, st%boxes(ib)%parts) ! ⚠ cutoff changes ?
            end if
          end block
        else
          block
            integer :: k
            call update_kappa_from_dist(st, 1)
            do k=1,1
              st%boxes(k)%ener = nr%energy_box(st, st%boxes(k))
            end do
          end block
        end if
      end if
    enddo
    call update_base_res(sim)
    rstep = step-first_step
    if(check_period(rstep, st%inp%step_log)) then
      call sgc_log_step(sim, step, output_unit)
    end if
    if(check_period(rstep, st%inp%step_nrj)) then
      block
        write(ukap,'(I10,ES16.7)') st%step, get_kappa(st,1)
        call output_psd_append(psd_interm, upsdi)
        write(upsdi, '(A)') "&" 
        ! Intermediate histogram -> reset at each output
        call histor8_reset(psd_interm)

        end block
    end if
    if(do_psd) then
      do i=1,st%boxes(1)%n
        call histor8_update(psd, st%boxes(1)%parts(i)%rayon)
      end do
    end if
    do i=1,st%boxes(1)%n
      call histor8_update(psd_interm, st%boxes(1)%parts(i)%rayon)
    end do
    call simulation_analysis_step(sim)
  end do
  call step_series_end(sim)
  close(ukap)
  close(upsdi)
  call write_move_statistics(sim,log_unit)
  write(log_unit, '(A15, F6.2,"% accepted, total trials: ", I0)') "Big steps: ", 100._8*real(big_accepted) / real(big_count), &
    big_count
  call simulation_analysis_finish(sim)
  call simulation_main_output(sim)
  call simulation_analysis_output(sim)
  call simulation_analysis_free(sim)
  if(log_unit /= output_unit) close(log_unit)
  end associate

end subroutine

real(8) function polydisp(st, ibox)
  type(PcmcState) :: st
  integer :: ibox
  integer :: i
  real(8) :: v,m
  v = 0.0
  m = 0.0
  do i=1,st%boxes(ibox)%n
    v = v + st%boxes(ibox)%parts(i)%rayon**2
    m = m + st%boxes(ibox)%parts(i)%rayon
  end do
  v = v / st%boxes(ibox)%n
  m = m / st%boxes(ibox)%n
  polydisp = (v-m*m)/m
end function

pure function mean_particle_volume(st, ibox) result(res)
  real(8) :: res
  type(PcmcState), intent(in) :: st
  integer, intent(in) :: ibox
  res = 4*PI/3 * sum(st%boxes(ibox)%parts(1:st%boxes(ibox)%n)%rayon**3)/st%boxes(ibox)%n
end function

pure function total_particle_volume(st, ibox) result(res)
  real(8) :: res
  type(PcmcState), intent(in) :: st
  integer, intent(in) :: ibox
  integer :: i
  res = 0.0_8
  do i = 1,st%boxes(ibox)%n
    res = res + st%boxes(ibox)%parts(i)%rayon**3
  end do
  res = res * (4._8*PI/3._8)
end function

function mean_particle_radius(st, ibox, r2m) result(res)
  real(8) :: res
  type(PcmcState), intent(in) :: st
  integer, intent(in) :: ibox
  real(8), intent(out) :: r2m
  integer :: i, np
  real(8) :: r
  np = 0
  res = 0.0_8
  r2m = 0.0_8
  do i = 1, st%boxes(ibox)%n
    r = st%boxes(ibox)%parts(i)%rayon
    ! if(r >= radius_min_toosmall) then
    !   np = np + 1
    !   res = res + r
    !   r2m = r2m + r*r
    ! end if
  end do
  if(np>0) then
    res = res / np
    r2m = r2m / np
  end if
end function


subroutine treat_additional_config(simu, sgcpars)
  type(MCsimulation) :: simu
  type(SgcParams) :: sgcpars
  character(len=:), allocatable :: val
  if(get_param(simu%state%config, "sgc.particle_resize_cvf_prob", val)) then
    read(val,*) sgcpars%particle_resize_cvf_prob
  end if
  if(get_param(simu%state%config, "sgc.particle_resize_prob", val)) then
    read(val,*) sgcpars%particle_resize_prob
  end if
  if(get_param(simu%state%config, "sgc.kappa_update_period", val)) then
    read(val,*) sgcpars%kappa_update_period
  end if
  if(get_param(simu%state%config, "sgc.big_step_kappa_update", val)) then
    sgcpars%big_mvmt = read_logical(val)
  end if
  if(get_param(simu%state%config, "sgc.max_radius", val)) then
    read(val,*) sgcpars%max_radius
  end if
  if(get_param(simu%state%config, "sgc.min_radius", val)) then
    read(val,*) sgcpars%min_radius
  end if
  if(get_param(simu%state%config, "sgc.max_radius_psd", val)) then
    read(val,*) sgcpars%max_radius_psd
  end if
  if(get_param(simu%state%config, "sgc.min_radius_psd", val)) then
    read(val,*) sgcpars%min_radius_psd
  end if
  if(get_param(simu%state%config, "sgc.resize_dr", val)) then
    read(val,*) sgcpars%dr_resize_max
  end if
  if(get_param(simu%state%config, "sgc.psd_dr", val)) then
    read(val,*) sgcpars%dr_psd
  end if
  if(get_param(simu%state%config, "sgc.npt_special", val)) then
    sgcpars%npt_special = read_logical(val)
  end if
  if(get_param(simu%state%config, "sgc.update_chempot_table", val)) then
    sgcpars%update_chempot_table = read_logical(val)
  end if
  if(get_param(simu%state%config, "sgc.chempot_update_factor", val)) then
    read(val,*) sgcpars%chempot_update_factor
  end if

  if(sgcpars%min_radius_psd > sgcpars%min_radius) sgcpars%min_radius_psd = sgcpars%min_radius
  if(sgcpars%max_radius_psd < sgcpars%max_radius) sgcpars%max_radius_psd = sgcpars%max_radius
end subroutine

integer function disappeared_particles(st, ib) result(res)
  type(PcmcState) :: st
  integer :: ib
  integer :: i
  res = 0
  do i=1,st%boxes(ib)%n
    if(st%boxes(ib)%parts(i)%rayon < sgcpars%min_radius) then
      res = res + 1
    end if
  end do
end function

integer function get_num_implicits(st, ib) result(res)
 type(PcmcState) :: st
  integer :: ib
  integer :: i
  res = 0
  do i=1,st%boxes(ib)%n
    if(st%boxes(ib)%parts(i)%numero < 0) then
      res = res + 1
    end if
  end do
end function

subroutine init_chempot_tables(ct, filename)
  use spline
  type(ChemPotData) :: ct
  character(*) :: filename
  integer :: f, n,i
  character(128) :: key
  real(8) :: rcell

  open(newunit=f, file=filename, status = 'old', action = 'read')

  read(f, *) key, n
  allocate(ct%radius(n), ct%table(n), ct%table_der2(n))
  do i = 1, n
    read(f, *) ct%radius(i), ct%table(i)
  end do
  
  do i=1, n
    call cubspline_init(n, ct%radius, ct%table, ct%table_der2) 
  end do
  close(f)
end subroutine

subroutine update_chempot(p, ct, distfile, mufile)
  type(HistoR8), intent(in) :: p
  type(ChemPotData), intent(in) :: ct
  character(*), intent(in) :: distfile, mufile

  integer :: i, fold, fmu, n
  real(8) :: dr, r, delta
  real(8), allocatable :: dist_goal(:), mu(:)
  real(8), parameter :: dmin = 1.e-30
  !call histor8_set_distribution(p)
  n = ubound(p%h, 1)
  if(n /= ubound(ct%radius, 1)) then
    write(error_unit, '(A)') "Error: Mean Distribution and μ table does not correspond, check min/max radii and bin"
    stop 1
  end if
  dr = 1._8 / p%dri
  do i = 1, n
    if( abs(ct%radius(i) - (i-0.5_8)*dr+p%rmin) < 1.d-6) then
      print *, "i = ", i, "r = ", ct%radius(i), "≠", (i-0.5_8)*dr+p%rmin
      write(error_unit, '(A)') "Error: Mean Distribution and μ table does not correspond, check min/max radii and bin"
      stop 1
    end if
  end do
  allocate(dist_goal(n), mu(n))

  open(newunit=fold, file=distfile, status="old")
  do i=1,n
    read(fold, *) r, dist_goal(i)
    if(abs(ct%radius(i) - (i-0.5_8)*dr+p%rmin) < 1.d-6) then
      write(error_unit, '(A,A,A)') "Error: Mismatch between μ table and goal distribution '",  distfile, "'"
      close(fold)
      stop 1
    end if
  end do
  close(fold)
  dist_goal = dist_goal / (sum(dist_goal)/dr)

  open(newunit=fold, file=mufile // ".old")
  write(fold, '(A, I0)') "#mu ", n
  do i=1,n
    write(fold, '(2E16.7)') ct%radius(i), ct%table(i) 
  end do
  close(fold)

  ! See Wilding 2009
  do i=1, n
    if(p%h(i) < dmin) then
      mu(i)= ct%table(i) + sgcpars%chempot_update_factor * log(dist_goal(i))/log(dmin)
    else
      !todo: use densities as volume is not always constant
      mu(i) = ct%table(i) + sgcpars%chempot_update_factor * log(dist_goal(i)/p%h(i))
    end if
  end do
  ! Shift (only Δμ is relevant for semigrand)
  mu(:) = mu(:) - mu(1)
  open(newunit=fmu, file=mufile)
  write(fmu, '(A, I0)') "#mu ", n
  do i=1,n
    write(fmu, '(2E16.7)') ct%radius(i), mu(i)
  end do
  print *, "TODEL ", sqrt(sum((mu(:) - ct%table(:))**2)/n)

  delta = sum(abs(p%h(:)-dist_goal(:)))
  write(output_unit, '(A, ES15.7)') "μ update: Δ = ", delta

  close(fmu)
end subroutine

subroutine init_tables(filename)
  use spline
  character(*) :: filename
  integer :: f, n1,n2, i,j
  character(128) :: key
  real(8) :: rcell

  open(newunit=f, file=filename, status = 'old', action = 'read')
  ! kappas
  read(f, *) key, n1
  call check(key, "kappa")
  allocate(kappa_table(n1))
  read(f, *) kappa_table
  ! radii
  read(f, *) key, n1
  call check(key, "radius")
  allocate(radius_table(n1))
  read(f, *) radius_table
  radius_table_min = minval(radius_table)
  radius_table_max = maxval(radius_table)
  ! Rcell
  read(f, *) key, n1, n2
  call check(key, "Rcell")
  allocate(rcell3_table(n1,n2))
  do j=1,n2
    do i=1,n1
      read(f,*) rcell
      rcell3_table(i,j) = rcell**3
    end do
  end do
  ! coeffs_charge
  read(f, *) key, n1, n2
  call check(key, "coeffs_charge")
  if(n1 /= sigma_kappa_n .or. n2 /= sigma_rad_n) then
    write(error_unit, *) "Error file '", filename, "': dimensions of 'coeffs_charge': got (", n1,&
      ",",n2,"), expected (",sigma_kappa_n, ",",sigma_rad_n,")"
    stop 1
  end if
  do j=1,n2
    do i=1,n1
      read(f,*) sigma_coeffs(i,j)
    end do
  end do
  close(f)

  ! Interpolate for each kappa
  block
    integer :: nr, nk
    nr = ubound(rcell3_table, 1)
    nk = ubound(rcell3_table, 2)
    allocate(rcell3_der2_table(nr, nk))
    do i=1, nk
      call cubspline_init(nr, radius_table, rcell3_table(:,i), rcell3_der2_table(:,i))
    end do

    ! check. To be removed later.
    do i=1, nk
      do j=1,nr
        if(abs(rcell3_table(j,i) - cubspline_eval(nr, radius_table, rcell3_table(:,i), &
          rcell3_der2_table(:,i), radius_table(j))) > 1.e-10*abs(rcell3_table(j,i))) then
          print *, rcell3_table(j,i), cubspline_eval(nr, radius_table, rcell3_table(:,i), &
          rcell3_der2_table(:,i), radius_table(j))
          stop "init_tables"
        end if
      end do
    end do
  end block

  contains

    subroutine check(a,b)
      use iso_fortran_env
      character(*) :: a,b
      if(trim(a) /= b) then
        write(error_unit, *) "Error file '", filename, "': Expected '",b, "', got '", trim(a), "'"
        stop 1
      end if
    end subroutine

end subroutine

real(8) function evalSurfaceCHarge(xi, kp) result(t)
  real(8), intent(in) :: xi, kp
  real(8) :: ai
  integer :: i,k
  t = 0._8
  do i=3,1,-1
    ai = sigma_coeffs(sigma_kappa_n, i)
    do k=sigma_kappa_n-1,1,-1
      ai = sigma_coeffs(k, i)  + kp*ai
    end do
    t =  ai + xi*t
  end do
end function 

real(8) function echarge(radius, kappa)
  !! Calc effective charge multiplied by e(kr)/(1+kr)
  real(8), intent(in) :: radius, kappa
  real(8) :: Z, xi

  if(radius <= radius_table_min) then
    xi = 1._8 / radius_table_min
  else if(radius >= radius_table_max) then
    xi = 1._8 / radius_table_max
  else
    xi = 1._8 / radius
  end if

  Z = 4* PI * radius**2 * evalSurfaceCHarge(xi, kappa)

  echarge = Z * exp(kappa * radius)/(1 + kappa*radius)
end function



subroutine update_kappa_from_dist(st, ibox)
  use spline
  !! Change the value of effective kappa for box `ibox`.
  type(PcmcState) :: st
  integer, intent(in) :: ibox

  integer :: i, p, nr, nk
  real(8), allocatable, save :: vcells(:), kappa2(:)
  real(8) :: new_kappa

  nk = ubound(kappa_table, 1)
  nr = ubound(radius_table, 1)

  if(.not. allocated(vcells)) then
    allocate(vcells(nk), kappa2(nk))
  end if

  Vcells = 0._8
  do p=1,st%boxes(ibox)%n
    do i=1,nk
      Vcells(i) = Vcells(i) + cubspline_eval(nr, radius_table, rcell3_table(:,i), &
        rcell3_der2_table(:,i), st%boxes(ibox)%parts(p)%rayon)
    end do
  end do
  Vcells =  - 4*PI/3 * Vcells ! Hack, because x values must be sorted by ascending order for spline

  
  call cubspline_init(nk, vcells, kappa_table, kappa2)
  if(- st%boxes(ibox)%volume < vcells(1)) then
    new_kappa = kappa_table(1)
  else
    new_kappa = cubspline_eval(nk, vcells, kappa_table, kappa2, - st%boxes(ibox)%volume)
  end if

  call set_kappa(st, ibox, new_kappa)

end subroutine

subroutine set_kappa(st, ibox, new_kappa)
  type(PcmcState) :: st
  integer, intent(in) :: ibox
  real(8), intent(in) :: new_kappa
  integer :: i

  ! Should it be here ?
  block
    real(8) :: old_kappa
    old_kappa = st%ndat%inv_debye_length
    st%ndat%kappa2 = new_kappa*new_kappa
    st%ndat%rcutoff_max = st%ndat%rcutoff_max * (old_kappa/new_kappa)
  end block

  do i=1,st%boxes(ibox)%n
    call change_part_echarge(st, ibox, i, echarge(st%boxes(ibox)%parts(i)%rayon, new_kappa))
  end do

  st%ndat%inv_debye_length = new_kappa

  do i=1,st%dist%nfam
    st%ndat%eff_charge(i) = echarge(st%dist%rad(i), new_kappa)
  end do
end subroutine

real(8) function get_kappa(st, ibox)
  type(PcmcState), intent(in) :: st
  integer, intent(in) :: ibox  
  get_kappa = st%ndat%inv_debye_length
end function

subroutine update_part_charge(st, ibox, i)
  type(PcmcState) :: st
  integer :: ibox !! box index
  integer :: i

  real(8) :: kappa
  kappa = st%ndat%inv_debye_length
  call change_part_echarge(st, ibox, i, echarge(st%boxes(ibox)%parts(i)%rayon, kappa))
end subroutine

subroutine  init_allcharges(st, ibox)
  !! Because convenience etc...
  type(PcmcState) :: st 
  integer :: ibox !! box index
  integer :: i
  real(8) :: kappa
  kappa = st%ndat%inv_debye_length
  do i=1, st%boxes(ibox)%n
    call change_part_echarge(st, ibox,i, echarge(st%boxes(ibox)%parts(i)%rayon, kappa))
  end do
end subroutine  init_allcharges

real(8) function chemPotent(ct, r) result(ret)
  use spline
  type(ChemPotData) :: ct
  real(8) :: r
  ! real(8), parameter ::  sigma = 1.2, rm = 8. ! sigma = 0.9625, rm = 13.75

  ! ret = 0._8
  ! ret = - 0.5*((r - rm)/sigma)**2 ! prefactor ?

  ret = cubspline_eval(ubound(ct%radius, 1), ct%radius, ct%table, ct%table_der2, r)
end function

integer function particle_resize(st, nr) result(ret)
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  integer :: ibox, i1
  real(8) :: r1p, r1pold, q1old, dbmu
  type(NRGRec) :: eold, de

  ret= mv_acc
  
  ibox = 1 ! ⚠ for now

  i1 = floor(ran2(st%rng)*st%boxes(ibox)%n) + 1

  r1pold = st%boxes(ibox)%parts(i1)%rayon
  q1old = st%boxes(ibox)%parts(i1)%ech

  r1p = r1pold + (ran2(st%rng)-0.5_8)*sgcpars%dr_resize_max
  if(r1p <= sgcpars%min_radius .or. r1p >= sgcpars%max_radius) then
    ret = mv_cir_rej
    return
  end if
  dbmu = chemPotent(chemdat, r1p) - chemPotent(chemdat, r1pold)   ! on oublie N(ro)/N(rn) et Λ(ro)/Λ(rn) 
  
  if(sgcpars%update_kappa_radius_resize) then
    eold = nr%energy_box(st,st%boxes(ibox))
    call change_part_radius(st, ibox, i1, r1p)
    call update_kappa_from_dist(st, ibox)
    de = nr%energy_box(st,st%boxes(ibox)) - eold
  else
    eold = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(i1))
    call change_part_radius(st, ibox, i1, r1p)
    call update_part_charge(st, ibox, i1)
    de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(i1)) - eold 
  end if

  if(de%noverlap) then
    if(metropolis(de%ep/kT - dbmu, st%rng)) then
      ret = mv_acc
      st%boxes(ibox)%ener%ep = st%boxes(ibox)%ener%ep + de%ep
    else
      ret = mv_nrj_rej
      call restore_state()
    end if
  else
    ret = mv_hs_rej
    call restore_state()
  end if

  contains

  subroutine restore_state()
    call change_part_radius(st, ibox, i1, r1pold)
    if(sgcpars%update_kappa_radius_resize) then
      call update_kappa_from_dist(st, ibox)
    else
      call change_part_echarge(st, ibox, i1, q1old)
    end if
  end subroutine

end function

integer function particle_resize_cvf(st, nr) result(ret)
  use iso_fortran_env
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  integer :: ibox, i1,i2
  real(8) :: r1p, r2p, r1pold, r2pold, q1old, q2old, dbmu
  type(NRGRec) :: eold1, eold2, de
  real(8), parameter :: onethrd = 1._8/3._8
  ret= mv_acc
  
  ibox = 1 ! ⚠ for now

  if(st%boxes(ibox)%n <= 1) then
    ret = mv_cir_rej
    return
  end if

  i1 = floor(ran2(st%rng)*st%boxes(ibox)%n) + 1

  i2 = floor(ran2(st%rng)*st%boxes(ibox)%n) + 1
  do while(i1 == i2)
    i2 = floor(ran2(st%rng)*st%boxes(ibox)%n) + 1
  end do

  r1pold = st%boxes(ibox)%parts(i1)%rayon
  r2pold = st%boxes(ibox)%parts(i2)%rayon
  q1old = st%boxes(ibox)%parts(i1)%ech
  q2old = st%boxes(ibox)%parts(i2)%ech

  r1p = r1pold + (ran2(st%rng)-0.5_8)*sgcpars%dr_resize_max
  if(r1p <= sgcpars%min_radius .or. r1p >= sgcpars%max_radius) then
    ret = mv_cir_rej
    return
  end if

  r2p = (r1pold**3 + r2pold**3 - r1p**3)
  if(r2p <= sgcpars%min_radius**3 .or. r2p > sgcpars%max_radius**3) then
    ret = mv_cir_rej
    return
  end if
    r2p = r2p**(onethrd)
  if(r2p == 0._8 .or. r1p == 0._8) then
    ! Particle 'disappear' ...
  end if
  
  dbmu = chemPotent(chemdat, r1p) + chemPotent(chemdat, r2p) - chemPotent(chemdat, r1pold) - chemPotent(chemdat, r2pold)
  
  ! Decom: ΔE1
  eold1 = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(i1))
  call change_part_radius(st, ibox, i1, r1p)
  call update_part_charge(st, ibox, i1)
  de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(i1)) - eold1

  ! ΔE2
  eold2 = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(i2))
  call change_part_radius(st, ibox, i2, r2p)
  call update_part_charge(st, ibox, i2)
  if(de%noverlap) de = de + (nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(i2)) - eold2)

  if(de%noverlap) then
    if(metropolis(de%ep/kT- dbmu, st%rng)) then
      ret = mv_acc
      st%boxes(ibox)%ener%ep = st%boxes(ibox)%ener%ep + de%ep
    else
      ret = mv_nrj_rej
    end if
  else
    ret = mv_hs_rej
  end if
  if(ret /= mv_acc) then
    call change_part_radius(st, ibox, i1, r1pold)
    call change_part_radius(st, ibox, i2, r2pold)
    call change_part_echarge(st, ibox, i1, q1old)
    call change_part_echarge(st, ibox, i2, q2old)
  end if

contains

end function

integer function changevol_npt_sgc(st, nr) result(ret)
  !! Volume change in the isobaric ensemble.
  implicit none
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  type(NrgRec) :: ener1_1,ener0_1
  integer :: ibox, truen
  real(8) :: deltabmu, Pref, oldkp
  real(8) :: fac1,dvol,rcell1_0(3),dumpv(3),oldbox(3,3)


  Pref = st%inp%Pref / (unitPsT*st%inp%temp) ! convert to kT units

  ibox = 1  ! ⚠ + floor(ran2(st%rng)*st%ntotbox)
  

  associate(labox => st%boxes(ibox))
    truen = labox%n

    ener0_1 = labox%ener

    dvol = (ran2(st%rng)-0.5_8)*st%inp%dv_vol
    if(st%inp%std_gibbs_volex) then
       fac1=1._8 + dvol/labox%volume
    else
       fac1 = exp(dvol)
    end if
    if(fac1<=0._8) then    
       ret=mv_cir_rej
       return
    end if
    if(st%inp%std_gibbs_volex) then
       deltabmu = - (truen)*log(fac1)  
    else
       deltabmu = - (truen+1)*log(fac1)
    end if
    deltabmu = deltabmu + Pref*labox%volume*(fac1-1._8)
    oldkp = get_kappa(st, ibox)
    oldbox = labox%tens
    call box_set_tensor(labox, fac1**(1._8/3._8)*labox%tens)
    call update_kappa_from_dist(st, ibox)
    if(st%celllist) then
       rcell1_0 = st%clistdecomp(ibox)%rcell
       dumpv(1)= st%rc/(labox%tens(1,1))
       dumpv(2)= st%rc/(labox%tens(2,2))
       dumpv(3)= st%rc/(labox%tens(3,3))
       call celldec_cree(st%clistdecomp(ibox),labox%n,labox%parts,&
            dumpv)
    endif
    ener1_1 = nr%energy_box(st,labox)
    if(ener1_1%noverlap) then
       if(metropolis((ener1_1%ep-ener0_1%ep)/kT + deltabmu,st%rng)) then
          labox%ener = ener1_1
          ret=mv_acc
       else
          ret=mv_nrj_rej
       end if
    else
      ret=mv_hs_rej
    end if
    if(ret /= mv_acc) then
      call box_set_tensor(labox, oldbox)
      call set_kappa(st, ibox, oldkp)
      if(st%celllist) then
        call celldec_cree(st%clistdecomp(ibox),labox%n,labox%parts,rcell1_0)
      end if
    end if
  end associate

end function

logical function check_energy(st, nr, ibox)
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  integer :: ibox
  type(NRGRec) :: e
  e = nr%energy_box(st, st%boxes(ibox))
  if(abs(st%boxes(ibox)%ener%ep-e%ep) >  abs(1.e-6*e%ep)) then
    check_energy = .false.
  else
    check_energy = .true.    
  end if
end function

subroutine output_psd(p, file, meanN)
  type(HistoR8) :: p
  real(8), intent(in) :: meanN
  character(*) :: file

  integer :: i, fo
  real(8) :: dr

  call histor8_set_distribution(p)
  open(newunit=fo, file=file)

  dr = 1._8 / p%dri
  do i=lbound(p%h,1),ubound(p%h,1)
    write(fo, '(2E16.7)') p%rmin + (i-0.5_8)*dr, p%h(i)*meanN
  end do

end subroutine

subroutine output_psd_append(p, unit)
  type(HistoR8) :: p
  integer :: unit

  integer :: i
  real(8) :: dr

  call histor8_set_distribution(p)

  dr = 1._8 / p%dri
  do i=lbound(p%h,1),ubound(p%h,1)
    write(unit, '(2E16.7)') p%rmin + (i-0.5_8)*dr, p%h(i)
  end do

end subroutine

end module
  
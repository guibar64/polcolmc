module ostwald
  !! Ostwald ripening with exchange with a finite bulk.
  
  !  Presently this module contains many code common with sgc, and a lot of cruft.
  use ieee_arithmetic
  use simulation
  use readparam
  use histograms
  use sgc, only: SgcParams
  implicit none

  type SolubData
    real(8) :: pKs = 2.7
    real(8) :: solid_density = 22.96 ! nm^-3
    real(8) :: lnKs

    ! surface tension
    real(8) :: gamma0 = 17 !(k_B*300K/nm^2)
    real(8) :: tolman_length = 1.0 ! nm
    real(8) :: gamma_doublelayer_inf = - 0.633177 
    real(8) :: gamma_doublelayer_slope = - 0.9722 / (- 0.633177)

    real(8) :: C = 22.96 ! prefactor (nm^-3)
    real(8) :: rmin = 1.0 ! nm

    real(8) :: kinetic_constant = 1.E-3      ! Unit may depend on way of
  end type

  type SolutionState
    integer(8) :: monos ! could be real, maybe
    real(8) :: vol !
  
    !integer(8) :: total_monos
  end type

  real(8), private ::particle_bulk_exchange_prob = 0.0, particle_bulk_formation_prob = 0.0
  real(8), private ::  min_radius = 2., max_radius = 30. ! not actually useful
  real(8), private ::  min_radius_psd = 2., max_radius_psd = 30., dr_psd = 0.10_8
  real(8), private ::  min_psd_rel = 0.05, max_psd_rel = 3.0, dr_psd_rel = 0.01_8
  integer, private :: steps_big_mvmt = 100
  logical, private :: big_mvmt = .true.
  logical, private :: npt_count_disappeared = .false.
  logical, private :: pseudo_kinetics = .false.


  real(8), allocatable :: kappa_table(:)  ! effective kappas
  real(8), allocatable :: radius_table(:) ! radii for R_cell interpolation
  real(8) :: radius_table_min, radius_table_max
  real(8), allocatable :: rcell3_table(:,:) ! R_cell(R, kappa)
  real(8), allocatable :: rcell3_der2_table(:,:)
  integer, parameter :: sigma_kappa_n = 4, sigma_rad_n = 3
  real(8) :: sigma_coeffs(sigma_kappa_n,sigma_rad_n) ! coeffs of polynomial fit of eff charge

  !real(8), allocatable :: kappacalc_vcells(:), kappacalc_kappa2(:)


  type(HistoR8), private :: psd, psd_interm, psd_renorm
  logical, private :: do_psd=.false.

  real(8), private :: dr_fam = 1.0_8, min_radius_fam = 6.0_8, max_radius_fam = 14.0_8
  integer, allocatable, private :: current_fams(:)
  integer, private :: gcmc_ntypes = -1

  integer(8), private :: minimum_monos_exchange=10, maximum_monos_exchange=100
  integer(8), private :: minimum_monos_particles = 10, maximum_monos_particles = 100
  logical, private :: bulk_exchange_create_destroy = .false.

  ! Solution details
  integer(8), private :: initial_number_dissolved_monos = -1_8
  integer(8), private :: total_number_monos = 300000_8
  logical, private :: finite_monomer_reservoir = .true.

  integer, private :: implicits_minimum_size = 2, implicits_maximum_size = 100 ! 100 ~ 1 nm for S
  real(8), private :: implicits_maximum_radius = 1.0
  real(8), private :: min_radius_gamma = 0.0 ! below that, γ(R) is constant
  integer, private :: partId = 0

  type(SolutionState), private :: solution
  type(SolubData), private :: solubd
  integer(8), allocatable, private :: number_monos(:) ! number of monomers for each particle

  !integer, allocatable :: implicits_num(:)   ! very small particles #parts(size)
  real(8), private :: nmonos_sol_mean, nmonos_count, nmonos_part_mean, nmonos_count_p

  integer(8), private :: part_created_count = 0_8, part_deleted_count = 0_8
  integer(8), private :: part_created_count_integ = 0_8, part_deleted_count_integ = 0_8

  real(8), allocatable, private :: deltaE_dissol(:)     ! Interaction part of the dissolution/precipitation ΔU
  real(8), allocatable, private :: deltaForm_dissol(:)  ! Internal part of the dissolution/precipitation ΔU
  integer(8), allocatable, private :: deltaU_dissol_count(:)
  real(8), allocatable, private :: deltaE_precip(:)     ! Interaction part of the dissolution/precipitation ΔU
  real(8), allocatable, private :: deltaForm_precip(:)  ! Internal part of the dissolution/precipitation ΔU
  integer(8), allocatable, private :: deltaU_precip_count(:)

  private :: nmonos_r8_to_rad, nmonos_to_rad, rad_to_nmonos, rdep_gamma, surface_nrg_raw, surface_nrg_r8, surface_nrg,&
    deltag_growth, deltag_part, deltamu, deltamu2, lnprefsol, ideal_aggregate_distribution, ideal_aggregate_distribution_n,&
    totally_not_kinetic_factor, calc_critical_radius, polydisp, mean_particle_radius,&
    disappeared_particles, get_num_implicits, evalSurfaceCharge, echarge, get_kappa, &
    vsgc_make_part, find_family_from_radius, particle_bulk_grow, particle_bulk_shrink, particle_mono_exchange, &
    particle_mono_exchange_vkin, particle_bulk_dissolution, particle_bulk_precipitation, particle_bulk_reformation, &
    changevol_npt_sgc, check_energy
  private :: initialize_solub, vsgc_log_step, run_simulation_vsgc, deltaU_precipdissol_reset, deltaU_precipdissol_print,&
    treat_additional_config, init_tables, update_kappa_from_dist, set_kappa, update_part_charge,&
    sac_add_part_box, sac_del_part_box, sac_change_part_size, output_psd, output_psd_append

contains

subroutine initialize_solub(sol, sd, monos, vol, rmin)
  type(SolutionState) :: sol
  type(SolubData) :: sd
  integer(8), intent(in) :: monos
  real(8), intent(in) :: vol, rmin
  sd%lnKs = -log(10.0_8)*sd%pKs + log(0.60221408_8) ! 1M/L = 0.60221408 nm-3
  sd%rmin = rmin
  sol = SolutionState(monos=monos, vol=vol)

  ! 🤔 Continuity. Well, there is no real model for small clusters, so...
  ! sd%C = exp(sd%lnKs + surface_nrg(sd, 1_8))
  sd%C = exp(sd%lnKs)
end subroutine

real(8) elemental function nmonos_r8_to_rad(sd, ns) result(ret)
  type(SolubData), intent(in) :: sd
  real(8), intent(in) :: ns
  ret = ((3._8*ns)/(4._8*PI*sd%solid_density))**(1._8/3._8)
end function

real(8) elemental function nmonos_to_rad(sd, ns) result(ret)
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: ns
  ret = nmonos_r8_to_rad(sd, real(ns, 8))
end function

integer(8) function rad_to_nmonos(sd, r) result(ret)
  type(SolubData) :: sd
  real(8) :: r
  ret = nint((4._8*PI*sd%solid_density*r**3)/3._8, 8)
end function

real(8) function rdep_gamma(sd, r) result(res)
  type(SolubData), intent(in) :: sd
  real(8), intent(in) :: r
  real(8) :: ri
  
  ri = 1/r
  res = (sd%gamma0 /(1+2*sd%tolman_length*ri) +&
      sd%gamma_doublelayer_inf*(1+sd%gamma_doublelayer_slope*ri))
end function

real(8) function surface_nrg_raw(sd, r) result(res)
  type(SolubData), intent(in) :: sd
  real(8), intent(in) :: r

  res =  4*PI*r*r*rdep_gamma(sd, r)
end function

real(8) function surface_nrg_r8(sd, r) result(res)
  type(SolubData), intent(in) :: sd
  real(8), intent(in) :: r
  res = surface_nrg_raw(sd, r)
  if(res <= 0.0 .or. r < sd%rmin) then
    res = 1.E7
  endif
end function

real(8) function surface_nrg(sd, ns) result(res)
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: ns

  real(8) :: r
  r = nmonos_to_rad(sd, ns)
  res = surface_nrg_r8(sd, r)
end function

function deltag_growth(sd, nsn, nso) result(res)
  real(8) :: res
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: nsn, nso
  integer(8) :: dn
  dn = nsn - nso  
  res = real(dn, 8) * (sd%lnKs) + surface_nrg(sd, nsn) - surface_nrg(sd, nso)
end function

function deltag_part(sd, ns, vol) result(res)
  real(8) :: res
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: ns
  real(8), intent(in) :: vol

  res = real(ns, 8) * (sd%lnKs) + surface_nrg(sd, ns) &
    - log(sd%C*vol)   ! maybe more a prefactor
end function

function deltamu(sd, ns, ca) result(res)
  real(8) :: res
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: ns
  real(8), intent(in) :: ca
  res = deltag_growth(sd, ns+1,ns) - log(ca)
end function

function deltamu2(sd, ns, na, v) result(res)
  real(8) :: res
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: ns, na
  real(8), intent(in) :: v
  res = deltag_growth(sd, ns+1,ns) - lnprefsol(-1_8, na, v)
end function

function lnprefsol(dn, n1at, vol) result(res)
  real(8) :: res
  ! type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: dn, n1at
  real(8), intent(in) :: vol
  real(8) :: pref, dnr, nar
  integer :: i
  pref = 1.0
  if(abs(dn) < 10_8) then
    if(dn > 0_8) then
      do i  = 1,dn
        pref = pref * vol / real(n1at + i, 8)
      end do
    else
      do i = 0, -dn-1
        pref = pref * real(n1at - i, 8) / vol
      end do
    end if
    res = log(pref)
  else
    dnr = real(dn, 8)
    nar = real(n1at, 8)
    res = dnr*log(vol) - log_gamma(nar+dnr) + log_gamma(nar)
  end if
end function

function ideal_aggregate_distribution(sd, r, c_a) result(res)
  real(8) :: res
  type(SolubData), intent(in) :: sd
  real(8), intent(in) :: r, c_a
  real(8) :: n
  n = real(rad_to_nmonos(sd, r), 8)
  res = 4*PI*sd%C*sd%solid_density* r*r*exp(n*(log(c_a) - sd%lnKs) - surface_nrg_r8(sd, r))
end function

function ideal_aggregate_distribution_n(sd, n, c_a) result(res)
  real(8) :: res
  type(SolubData), intent(in) :: sd
  real(8), intent(in) :: n, c_a
  real(8) :: r
  r =  nmonos_r8_to_rad(sd, n)
  res = exp(n*(log(c_a) - sd%lnKs) - surface_nrg_r8(sd, r) &
   !+ sd%lnKs + surface_nrg_r8(sd, r))
   )*sd%C
end function

real(8) function totally_not_kinetic_factor(sd, n) result(ret)
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: n
  real(8) :: r
  r = nmonos_to_rad(sd, n)
  ! alpha = sd%gamma0/sd%solid_density
  ret = sd%kinetic_constant * r**2 
end function

real(8) function calc_critical_radius(sd, n_min, n_max, c_monos) result(ret)
  type(SolubData), intent(in) :: sd
  integer(8), intent(in) :: n_min, n_max
  real(8), intent(in) :: c_monos
  integer(8) :: low, high, n
  low = n_min
  high = n_max
  do while(high-low>1_8)
    n=(high+low) / 2_8
    ! ! Supposes Δμ ↓ with n
    if(deltamu(sd, n, c_monos) < 0.0_8) then
      high = n
    else
      low=n
    end if
  end do
  ret = 0.5_8*(nmonos_to_rad(sd, low) + nmonos_to_rad(sd, high))
end function


subroutine ostwald_main(simu)
  use iso_fortran_env
  type(MCsimulation) :: simu
  logical :: oldprod
  character(len=:), allocatable :: prefix
  real(8) :: prob_volmod

  if(simu%state%ntotbox > 1) then
    write(output_unit, '(A)') "This simulation is not compatible with multiple boxes. Exiting"
    stop 1
  end if 

  call treat_additional_config(simu)


  call add_move(simu%moves, "Particle growth (bulk)" , 0.5_8*particle_bulk_exchange_prob, particle_bulk_grow)
  call add_move(simu%moves, "Particle shrinking (bulk)" , 0.5_8*particle_bulk_exchange_prob, particle_bulk_shrink)

  call add_move(simu%moves, "Particle precipitation" , 0.5_8*particle_bulk_formation_prob, particle_bulk_precipitation)
  call add_move(simu%moves, "Particle dissolution" , 0.5_8*particle_bulk_formation_prob, particle_bulk_dissolution)

  if(simu%state%inp%asimu == "Ostwald-PT") then
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

  ! Init state distribution
  block
    integer :: nfam, j, ib = 1, i
    nfam = floor((max_radius_fam-min_radius_fam)/dr_fam) + 1
    gcmc_ntypes = nfam !(max_radius_fam-min_radius_fam)/dr_fam
    allocate(current_fams(nfam))
    current_fams = 0
    associate(b => simu%state%boxes(ib))
    do i=1, b%n
      j = find_family_from_radius(simu%state, b%parts(i)%rayon)
      current_fams(j) = current_fams(j) + 1
    end do
    end associate
  end block


  block
    integer :: i
    allocate(number_monos(ubound(simu%state%boxes(1)%parts,1)))

    !call print_clist(output_unit,simu%state%clistdecomp(1),simu%state%boxes(1)%n,simu%state%boxes(1)%parts)
    do i=1, simu%state%boxes(1)%n
      ! print *, "ahhh:: ", simu%state%boxes(1)%parts(i)%rayon
      number_monos(i) = rad_to_nmonos(solubd, simu%state%boxes(1)%parts(i)%rayon)
      ! Rounding error 🤔
      simu%state%boxes(1)%parts(i)%rayon = nmonos_to_rad(solubd, number_monos(i))

      if(number_monos(i) <= implicits_maximum_size) then
        ! simu%state%boxes(1)%parts(i)%numero = -1 ! out of the loop
        simu%state%boxes(1)%parts(i)%ech = 0.0_8
      end if
    end do
    partId = simu%state%boxes(1)%n + 1

  end block
  block

    ! volume -> excess volume
    write(*,'(A, G15.8)') "Total number of particle monomers = ", sum(number_monos(1:simu%state%boxes(1)%n))
    if(initial_number_dissolved_monos < 0_8) &
      initial_number_dissolved_monos = total_number_monos - sum(number_monos(1:simu%state%boxes(1)%n))
    call initialize_solub(solution, solubd, initial_number_dissolved_monos, simu%state%boxes(1)%volume, min_radius_gamma)

    ! Initialize Implicit small
    !call implicits_init(implicits_num, implicits_minimum_size, implicits_maximum_size)
    implicits_maximum_radius = nmonos_to_rad(solubd, int(implicits_maximum_size,8))
  end block
    ! Initialize size histogram
  block
    integer :: npsd, i, fps
    real(8) :: r, npid, psdid
    integer(8) :: j
    npsd = floor((max_radius_psd - min_radius_psd) / dr_psd) - 1
    call histor8_init(psd, dr_psd, min_radius_psd, npsd)
    do i=1, simu%state%boxes(1)%n
      call histor8_update(psd, simu%state%boxes(1)%parts(i)%rayon)
    end do
  

    call histor8_init(psd_interm, dr_psd, min_radius_psd, npsd)
    do i=1, simu%state%boxes(1)%n
      call histor8_update(psd_interm, simu%state%boxes(1)%parts(i)%rayon)
    end do

    npsd = floor((max_psd_rel - min_psd_rel) / dr_psd_rel) - 1
    call histor8_init(psd_renorm, dr_psd_rel, min_psd_rel, npsd)

    nmonos_count = 0.0_8
    nmonos_sol_mean = 0.0_8
    nmonos_count_p = 0.0_8
    nmonos_part_mean = 0.0_8


    open(newunit=fps, file="pcmc_gamma.dat")
    do i=1,npsd
      r = min_radius_psd+(i+0.5_8)*dr_psd
      psdid = rdep_gamma(solubd, r)
      write(fps, '(2ES16.7)') r, psdid
    end do
    close(fps)
  end block
  ! Comp. PSD for prod run only
  do_psd=.false.

  !stop "do not run simulation"
  oldprod = simu%prod
  simu%prod = .false.
  prefix = simu%prefix
  simu%prefix = prefix // "_eq"
  call run_simulation_vsgc(simu, simu%state%inp%maxcycle_stab)

  do_psd = .true.
  simu%prod = oldprod
  simu%prefix = prefix
  call run_simulation_vsgc(simu, simu%state%inp%maxcycle_calc)
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

    nmonos_sol_mean = nmonos_sol_mean / nmonos_count
    nmonos_part_mean = nmonos_part_mean / nmonos_count_p
    write(output_unit, '("<monos sol> = ",G15.7," <c_a> = ",G15.7," mol/L")') nmonos_sol_mean, &
      nmonos_sol_mean / solution%vol / 0.60221408
    write(output_unit, '("<monos part> = ",G15.7," (",G15.7," nm)")') nmonos_part_mean, &
      nmonos_r8_to_rad(solubd, nmonos_part_mean)
    ! print *, part_created_count / real(simu%state%inp%maxcycle_calc + simu%state%inp%maxcycle_stab)
  
    print '("Main: Created/Destroyed      ", I0,"/",I0)', part_created_count, part_deleted_count
    print '("Separate: Created/Destroyed  ", I0,"/",I0)', part_created_count_integ, part_deleted_count_integ
  end block
end subroutine

subroutine vsgc_log_step(sim, step, log_unit)
  type(MCSimulation) :: sim
  integer, intent(in) :: step, log_unit
  call log_step(sim, step, log_unit)
  associate(st => sim%state)

  write(log_unit,'(A,G15.7, A)') "kappa                        = ", get_kappa(st,1), " nm⁻¹"
  write(log_unit,'(A,I0)') "#implicits                = ", get_num_implicits(st, 1)
  ! write(log_unit, '(A,G15.7,A,I0, A,G15.7,A)') "V_p = ", total_particle_volume(st, 1), " = 4π/3 × ", st%boxes(1)%n, &
  !   " × ", (mean_particle_volume(st, 1)*3/(4*PI))**(1._8/3._8),"^3"
  write(log_unit, '(A,I0,A,I0)') "#monos = ", solution%monos + sum(number_monos(1:st%boxes(1)%n)), " sol = ", solution%monos
  write(log_unit, '(A, G15.7)') "Ca / Ks = ", solution%monos/st%boxes(1)%volume*exp(-solubd%lnKs)
  print '("Δμ = ", G16.7, " R_c = ", G16.7)', deltamu2(solubd, number_monos(1), solution%monos, st%boxes(1)%volume),&
    calc_critical_radius(solubd, 2_8, rad_to_nmonos(solubd, max_radius_psd), solution%monos/st%boxes(1)%volume)
      !deltamu2(sd, ns, na, v)

  end associate
end subroutine

subroutine run_simulation_vsgc(sim, maxsteps)
  !! Runs a simulation
  type(MCSimulation) :: sim
  integer, intent(in) :: maxsteps
  integer :: first_step, step, rstep, log_unit
  integer :: sel_mvmt, ins, inner_count, i, big_count, big_accepted, smallsteps_count
  integer :: ininnerc, ukap, urmean, upsdi, upsdrel
  real(8) :: r2mean
  type(Box) :: oldbox
  integer(8), allocatable :: old_number_monos(:)
  type(SolutionState) :: old_solution

  !open(newunit=log_unit, file=sim%prefix // "_log.txt")
  log_unit = output_unit
  associate(st => sim%state, nr => sim%nrg)
  call reset_all_moves(sim%moves) ! Maybe not necessary
  first_step = sim%state%step
  call init_step_series(sim,sim%prefix // "_traj.inst", sim%prefix,&
        first_step=first_step)
  call simulation_analysis_init(sim)
  call vsgc_log_step(sim, first_step, output_unit)
  
  open(newunit=ukap, file = sim%prefix // "_kappa.dat")
  open(newunit=urmean, file = sim%prefix // "_rmean.dat")
  block
    real(8) :: rmean, r2mean
    rmean = mean_particle_radius(st, 1, r2mean)
    write(urmean,'(I10,4ES16.7)') st%step, rmean, &
      calc_critical_radius(solubd, 2_8, rad_to_nmonos(solubd, max_radius_psd), solution%monos/solution%vol), &
      r2mean-rmean*rmean, real(st%step*st%sdata%nnn)/st%boxes(1)%n
  end block

  open(newunit=upsdi, file = sim%prefix // "_intermpsd.dat")
  open(newunit=upsdrel, file = sim%prefix // "_psdrel.dat")
  if(st%inp%internal_steps > 0) then
    st%sdata%nnn = st%inp%internal_steps
  else
    st%sdata%nnn = st%boxes(1)%n
  end if

  if(steps_big_mvmt > 0) then
    inner_count = sim%state%sdata%nnn / steps_big_mvmt
    smallsteps_count = steps_big_mvmt
  else
    inner_count = sim%state%sdata%nnn
    smallsteps_count = 1
  end if
  print *, "internal_steps = ", sim%state%sdata%nnn
  big_count = 0
  big_accepted = 0
  call deltaU_precipdissol_reset(sim%state)
  oldbox = st%boxes(1) !
  allocate(old_number_monos, source=number_monos)
  old_solution = solution
  do step = 1 + first_step, first_step + maxsteps
    sim%state%step = step ! 🤔
    call update_trans_ampl(sim%state)
    do ins=1, inner_count
      if(steps_big_mvmt > 0 .and. ins == inner_count) then
        ininnerc = modulo(sim%state%sdata%nnn, steps_big_mvmt)
      else
        ininnerc = smallsteps_count
      end if
      do i = 1, smallsteps_count
        ininnerc = ininnerc + 1
        sel_mvmt = select_move(sim%moves, ran2(sim%state%rng))
        call new_movinst(sim%state, sim%nrg, sim%moves%mv(sel_mvmt))
      end do

      if(steps_big_mvmt > 0) then
        if(big_mvmt) then
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
            !epn = nr%energy_box(st, st%boxes(ib))
            epn =  st%boxes(ib)%ener
            call update_kappa_from_dist(st, ib) ! ⚠ cutoff changes ?
            en = nr%energy_box(st, st%boxes(ib))
            ! p_ref = st%inp%Pref / (unitPsT*st%inp%temp)
            ! if(st%inp%asimu(1:3) == "Ostwald-NPT") then
            !   dpv = p_ref*(st%boxes(ib)%volume-oldbox%volume)
            ! else
            !   dpv = 0.0_8
            ! end if
            if(en%noverlap .and. metropolis((en%ep-epn%ep)/kT, st%rng)) then
              big_accepted = big_accepted + 1
              st%boxes(ib)%ener = en
              oldbox = st%boxes(ib)
              old_number_monos = number_monos
              old_solution = solution
            else
              call set_kappa(st, ib, oldkp)
              st%boxes(ib) = oldbox
              number_monos = old_number_monos
              solution = old_solution
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
      call vsgc_log_step(sim, step, output_unit)
    end if
    if(check_period(rstep, st%inp%step_nrj)) then
      solution%vol = st%boxes(1)%volume ! should be updated after every volume change, but...
      block
        real(8) :: rmean, r2mean
        write(ukap,'(I10,ES16.7)') st%step, get_kappa(st,1)

        rmean = mean_particle_radius(st, 1, r2mean)
        write(urmean,'(I10,4ES16.7)') st%step, rmean, &
          calc_critical_radius(solubd, 2_8, rad_to_nmonos(solubd, max_radius_psd), solution%monos/solution%vol), &
          r2mean-rmean*rmean, real(st%step*st%sdata%nnn)/st%boxes(1)%n

        call output_psd_append(psd_interm, upsdi)
        write(upsdi, '(A)') "&" 
        ! Intermediate histogram -> reset at each output
        call histor8_reset(psd_interm)

        ! Instantaneous normalized psd
        call histor8_reset(psd_renorm)
        do i=1,st%boxes(1)%n
          call histor8_update(psd_renorm, st%boxes(1)%parts(i)%rayon/rmean)
        end do
        call output_psd_append(psd_renorm, upsdrel)
          write(upsdrel, '(A)') "&" 
        end block
    end if
    if(do_psd) then
      do i=1,st%boxes(1)%n
        call histor8_update(psd, st%boxes(1)%parts(i)%rayon)
      end do
      nmonos_count = nmonos_count + 1._8
      nmonos_sol_mean = nmonos_sol_mean + solution%monos
      nmonos_count_p = nmonos_count_p + st%boxes(1)%n
      nmonos_part_mean = nmonos_part_mean + sum(number_monos(1:st%boxes(1)%n))
    end if
    do i=1,st%boxes(1)%n
      call histor8_update(psd_interm, st%boxes(1)%parts(i)%rayon)
    end do
    call simulation_analysis_step(sim)
  end do
  call step_series_end(sim)
  close(ukap)
  close(urmean)
  close(upsdi)
  call write_move_statistics(sim,log_unit)
  write(log_unit, '(A15, F6.2,"% accepted, total trials: ", I0)') "Big steps: ", 100._8*real(big_accepted) / real(big_count), &
    big_count
  call deltaU_precipdissol_print(sim%state, log_unit, sim%prefix // "_dFormvsdEint.dat")
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

! pure function mean_particle_volume(st, ibox) result(res)
!   real(8) :: res
!   type(PcmcState), intent(in) :: st
!   integer, intent(in) :: ibox
!   res = 4*PI/3 * sum(st%boxes(ibox)%parts(1:st%boxes(ibox)%n)%rayon**3)/st%boxes(ibox)%n
! end function

! pure function total_particle_volume(st, ibox) result(res)
!   real(8) :: res
!   type(PcmcState), intent(in) :: st
!   integer, intent(in) :: ibox
!   integer :: i
!   res = 0.0_8
!   do i = 1,st%boxes(ibox)%n
!     res = res + st%boxes(ibox)%parts(i)%rayon**3
!   end do
!   res = res * (4._8*PI/3._8)
! end function

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
    if(number_monos(i) > implicits_maximum_size) then
      np = np + 1
      res = res + r
      r2m = r2m + r*r
    end if
  end do
  if(np>0) then
    res = res / np
    r2m = r2m / np
  end if
end function

! function total_particle_volume_monos(st, ibox) result(res)
!   real(8) :: res
!   type(PcmcState), intent(in) :: st
!   integer, intent(in) :: ibox
!   integer :: i
!   res = 0.0_8
!   do i = 1,st%boxes(ibox)%n
!     res = res + 4*PI/3 * nmonos_to_rad(solubd, number_monos(i))**3
!   end do
!   do i = 1,st%boxes(ibox)%n
!     res = res - 4*PI/3 * st%boxes(ibox)%parts(i)%rayon**3
!   end do
! end function

subroutine deltaU_precipdissol_reset(st)
  type(PcmcState) :: st
  integer :: n
  n = st%dist%nfam
  if(.not. allocated(deltaE_dissol)) then
    allocate(deltaE_dissol(n), deltaForm_dissol(n), deltaE_precip(n), deltaForm_precip(n), deltaU_dissol_count(n),&
      deltaU_precip_count(n))
  end if
  deltaE_dissol = 0._8   
  deltaForm_dissol = 0._8
  deltaU_dissol_count = 0
  deltaE_precip = 0._8    
  deltaForm_precip = 0._8
  deltaU_precip_count = 0
end subroutine


subroutine deltaU_precipdissol_print(st, fu, file)
  type(PcmcState), intent(in) :: st
  integer, intent(in) :: fu
  character(*), intent(in) :: file
  integer :: fo, f
  real(8) :: dep, ded, dfp, dfd
  write(fu, '("<ΔE>    dissol=", ES16.7," precip=",ES16.7)') sum(deltaE_dissol)/sum(deltaU_dissol_count), &
    sum(deltaE_precip)/sum(deltaU_precip_count)
  write(fu, '("<ΔForm> dissol=", ES16.7," precip=",ES16.7)') sum(deltaForm_dissol)/sum(deltaU_dissol_count), &
    sum(deltaForm_precip)/sum(deltaU_precip_count)

  open(newunit=fo, file=file)
  write(fo, '(A)') "# R ΔE_int_diss ΔF_form_diss ΔE_int_prec ΔF_form_prec"
  do f=1,st%dist%nfam
    if(deltaU_dissol_count(f) > 0) then
      ded = deltaE_dissol(f)/deltaU_dissol_count(f)
      dfd = deltaForm_dissol(f)/deltaU_dissol_count(f)
    else
      ded = 0._8
      dfd = 0._8
    end if
    if(deltaU_precip_count(f) > 0) then
      dep = deltaE_precip(f)/deltaU_precip_count(f)
      dfp = deltaForm_precip(f)/deltaU_precip_count(f)
    else
      dep = 0._8
      dfp = 0._8
    end if
    write(fo, '(5ES16.8)') st%dist%rad(f), ded, dfd, dep, dfp
  end do
  close(fo)
end subroutine

subroutine treat_additional_config(simu)
  type(MCsimulation) :: simu
  character(len=:), allocatable :: val

  if(get_param(simu%state%config, "ostwald.max_radius", val)) then
    read(val,*) max_radius
  end if
  if(get_param(simu%state%config, "ostwald.min_radius", val)) then
    read(val,*) min_radius
  end if
  if(get_param(simu%state%config, "ostwald.max_radius_psd", val)) then
    read(val,*) max_radius_psd
  end if
  if(get_param(simu%state%config, "ostwald.min_radius_psd", val)) then
    read(val,*) min_radius_psd
  end if
  if(get_param(simu%state%config, "ostwald.big_step_kappa_update", val)) then
    big_mvmt = read_logical(val)
  end if
  if(get_param(simu%state%config, "ostwald.kappa_update_period", val)) then
    read(val,*) steps_big_mvmt
  end if

  if(get_param(simu%state%config, "ostwald.particle_bulk_exchange_prob", val)) then
    read(val,*) particle_bulk_exchange_prob
  end if
  if(get_param(simu%state%config, "ostwald.bulk_exchange_create_destroy", val)) then
    read(val,*) bulk_exchange_create_destroy
  end if

  if(get_param(simu%state%config, "ostwald.minimum_monos_exchange", val)) then
    read(val,*) minimum_monos_exchange
  end if
  if(get_param(simu%state%config, "ostwald.maximum_monos_exchange", val)) then
    read(val,*) maximum_monos_exchange
  end if
  if(get_param(simu%state%config, "ostwald.initial_number_dissolved_monos", val)) then
    read(val,*) initial_number_dissolved_monos
  end if
  if(get_param(simu%state%config, "ostwald.total_number_monos", val)) then
    read(val,*) total_number_monos
  end if

  if(get_param(simu%state%config, "ostwald.pKs", val)) then
    read(val,*) solubd%pKs
  end if
  if(get_param(simu%state%config, "ostwald.solid_density", val)) then
    read(val,*) solubd%solid_density
  end if
  if(get_param(simu%state%config, "ostwald.gamma0", val)) then
    read(val,*) solubd%gamma0
  end if
  if(get_param(simu%state%config, "ostwald.gamma_doublelayer_inf", val)) then
    read(val,*) solubd%gamma_doublelayer_inf
  end if
  if(get_param(simu%state%config, "ostwald.gamma_doublelayer_slope", val)) then
    read(val,*) solubd%gamma_doublelayer_slope
  end if
  if(get_param(simu%state%config, "ostwald.tolman_length", val)) then
    read(val,*) solubd%tolman_length
  end if

  if(get_param(simu%state%config, "ostwald.particle_bulk_formation_prob", val)) then
    read(val,*) particle_bulk_formation_prob
  end if
  if(get_param(simu%state%config, "ostwald.finite_monomer_reservoir", val)) then
    read(val,*) finite_monomer_reservoir
  end if

  if(get_param(simu%state%config, "ostwald.implicits_minimum_size", val)) then
    read(val,*) implicits_minimum_size
  end if
  if(get_param(simu%state%config, "ostwald.implicits_maximum_size", val)) then
    read(val,*) implicits_maximum_size
  end if
  if(get_param(simu%state%config, "ostwald.pseudo_kinetics", val)) then
    read(val,*) pseudo_kinetics
  end if

  if(get_param(simu%state%config, "ostwald.minimum_monos_particles", val)) then
    read(val,*) minimum_monos_particles
  end if
  if(get_param(simu%state%config, "ostwald.min_radius_gamma", val)) then
    read(val,*) min_radius_gamma
  end if

  if(get_param(simu%state%config, "ostwald.maximum_monos_particles", val)) then
    read(val,*) maximum_monos_particles
  end if

  if(min_radius_psd > min_radius) min_radius_psd = min_radius
  if(max_radius_psd < max_radius) max_radius_psd = max_radius
end subroutine

integer function disappeared_particles(st, ib) result(res)
  type(PcmcState) :: st
  integer :: ib
  integer :: i
  res = 0
  do i=1,st%boxes(ib)%n
    if(st%boxes(ib)%parts(i)%rayon < min_radius) then
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

! real(8) function echarge(radius, kappa)
!   !! Calc effective charge multiplied by e(kr)/(1+kr)
!   real(8), intent(in) :: radius, kappa
!   real(8) :: Z

!   Z = 4* PI * radius**2 * sigmaPlan * (1+ Acharge/(kappa * radius))

!   echarge = Z * exp(kappa * radius)/(1+ kappa*radius)
! end function

! subroutine change_volume(st, ibox, factor)
!   !! Updates quantities affected by a volume change...
!   type(PcmcState) :: st
!   integer :: ibox !! box index
!   real(8) :: factor !! volume change factor

!   integer :: i
!   real(8) :: kappa, phi

!   phi = 0.0
!   do i=1,st%boxes(ibox)%n
!     phi = phi + st%boxes(ibox)%parts(i)%rayon**3
!   end do
!   phi = phi * 4.0*PI/3.0 / st%boxes(ibox)%volume
  
!   ! Just for testing...
!   kappa = kappa_ref *(phi/phi_ref)**2

!   !! ⚠ It affects all boxes !
!   st%ndat%inv_debye_length = kappa

!   do i=1,st%boxes(ibox)%n
!     st%boxes(ibox)%parts(i)%ech = echarge(st%boxes(ibox)%parts(i)%rayon, kappa)
!   end do

! end subroutine

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


function vsgc_make_part(st, pos, rad) result(part)
  type(Particule) :: part
  type(PcmcState), intent(in) :: st
  real(8), intent(in) :: pos(3)
  real(8), intent(in) :: rad
  part%pos = pos
  part%rayon = rad
  part%ech = echarge(rad, st%ndat%inv_debye_length)
  part%famille = find_family_from_radius(st, rad)
end function


integer function find_family_from_radius(st, r) result(fam)
  ! Updates family within the predefined distribution.
  ! ⚠ Supposes a sorted distribution, may not work.
  type(PcmcState) :: st
  real(8), intent(in) :: r
  integer :: j
  real(8) :: mr
  ! 🤔 Pb of interval vs center
  j = floor((r + 1.e-10_8 - min_radius_fam)/dr_fam) + 1
  if(j<1) then
    fam = 1
  else if(j > ubound(current_fams,1)) then
    fam = ubound(current_fams,1)
  else
    fam = j
  end if
end function

subroutine sac_add_part_box(st, ibox, part, n)
  type(PcmcState), intent(inout) :: st
  integer, intent(in) :: ibox
  type(Particule), intent(inout) :: part
  integer(8) :: n
  integer, save :: counter = 0
  integer :: np
  counter = counter + 1
  np = st%boxes(ibox)%n
  if(n<=implicits_maximum_size) then
    ! part%numero = -1 ! out of the loop, sort of
    part%ech = 0.0
  end if
    part%numero = partId ! 🤔 there will be overlap
    if(partId >= huge(1)) then
      partId = np+1
    else
      partId = partId + 1      
    end if
  ! end if
  part%rayon = nmonos_to_rad(solubd, n)
  call add_part_box(st, ibox, part)
  if(np+1 <= ubound(number_monos,1)) then
    number_monos(np+1) = n
  else
    stop "SAC ADD: Too many particles"
  end if

end subroutine

subroutine sac_del_part_box(st, ibox, ip)
  type(PcmcState), intent(inout) :: st
  integer, intent(in) :: ibox
  integer, intent(in) :: ip
  integer(8) :: n
  integer :: i, nid
  nid = st%boxes(ibox)%parts(ip)%numero

  call remove_part_box(st, ibox, ip)
  do i=ip, st%boxes(ibox)%n
    number_monos(i) = number_monos(i+1)
  end do

  
end subroutine

subroutine sac_change_part_size(st, ibox, ip, n)
  type(PcmcState) :: st
  integer, intent(in) :: ibox
  integer, intent(in) :: ip
  integer(8), intent(in) :: n
  real(8) :: r
  integer :: tip
  r = nmonos_to_rad(solubd, n)
  call change_part_radius(st, ibox, ip, r)
  call update_part_charge(st, ibox, ip)

  if(n<=implicits_maximum_size) then
    ! Desactivate particle
    st%boxes(ibox)%parts(ip)%ech = 0.0_8    
  else
    ! Activate particle
  end if
  number_monos(ip) = n

end subroutine


integer function particle_bulk_grow(st, nr) result(ret)
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  if(pseudo_kinetics) then
    ret = particle_mono_exchange_vkin(st, nr, dissolution=.false.)
  else
    ret = particle_mono_exchange(st, nr, dissolution=.false.)
  end if

end function
integer function particle_bulk_shrink(st, nr) result(ret)
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  if(pseudo_kinetics) then
    ret = particle_mono_exchange_vkin(st, nr, dissolution=.true.)
  else
    ret = particle_mono_exchange(st, nr, dissolution=.true.)
  end if
end function

integer function particle_mono_exchange(st, nr, dissolution) result(ret)
  use iso_fortran_env
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  logical, intent(in) :: dissolution
  integer :: ibox, p,k, nsmall, j, i, np_tot, f
  integer(8) :: dn, m0, m1, nso
  real(8) :: r1p, r2p, fac, r1pold, r2pold, q1old, q2old, q1p, &
    q2p, pos(3), vpo, vpn, lpref, dg
  type(NRGRec) :: eold, eold2, de
  real(8), parameter :: onethrd = 1._8/3._8
  type(Particule) :: part
  logical :: addordelpart
  integer, save :: counter = 0
  counter = counter + 1
  ret= mv_acc
  ibox = 1 ! ⚠ for now
  associate(lbox => st%boxes(ibox))
  solution%vol = lbox%volume ! should be updated after every volume change, but...
  addordelpart = .false.
  dn = minimum_monos_exchange + floor(ran2(st%rng) * (maximum_monos_exchange-minimum_monos_exchange))
  np_tot = st%boxes(ibox)%n
  ! vpo = total_particle_volume(st, ibox)
  if(dissolution) then
    ! Dissolution
    if(np_tot > 0) then
      ! selects a particule
      p = floor(ran2(st%rng)*np_tot) + 1
      if(number_monos(p) >= dn + minimum_monos_exchange) then
        f = st%boxes(ibox)%parts(p)%famille
        dg = deltag_growth(solubd, number_monos(p)-dn, number_monos(p))
        lpref = lnprefsol(dn, solution%monos, solution%vol)
        if(bulk_exchange_create_destroy) lpref = lpref + log(real(lbox%n,8)/real(lbox%n+1, 8))
        !lpref = lpref + log(totally_not_kinetic_factor(solubd, number_monos(p)))

        nso = number_monos(p)
        eold = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        call sac_change_part_size(st, ibox, p, nso-dn)
        de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p)) - eold

        ret = acceptsornot(de, lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = dn + solution%monos
        else
          call sac_change_part_size(st, ibox, p, nso)
        end if
        ! only if acc ?
        deltaForm_dissol(f) = deltaForm_dissol(f) + (dg-lpref)/dn
        deltaE_dissol(f) = deltaE_dissol(f) + de%ep/dn
        deltaU_dissol_count(f) = deltaU_dissol_count(f) + 1

      else if(bulk_exchange_create_destroy .and. number_monos(p) == dn) then
        ! delete the particle
        addordelpart = .true.
        dg = - deltag_part(solubd, dn, lbox%volume)
        lpref = lnprefsol(dn, solution%monos, solution%vol)
        part = lbox%parts(p)
        de = - nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        call sac_del_part_box(st, ibox, p)

        ret = acceptsornot(de, lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = dn + solution%monos
          part_deleted_count = part_deleted_count + 1
        else
          call sac_add_part_box(st, ibox, part, dn)
        end if
      else
        ret = mv_cir_rej
      end if
    else
      ret = mv_cir_rej
    end if
  else
    ! Precipitation
    if(solution%monos >= dn) then
      if(bulk_exchange_create_destroy) then
        p = floor(ran2(st%rng)*(np_tot+1)) + 1
      else
        p = floor(ran2(st%rng)*(np_tot)) + 1
      end if
      ! p = floor(ran2(st%rng)*(st%boxes(ibox)%n)) + 1
      if(p <= np_tot) then
        f = st%boxes(ibox)%parts(p)%famille
        nso = number_monos(p)
        !call test_over(st, "Precip: Grow Before")
        eold = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        call sac_change_part_size(st, ibox, p, nso+dn)
        de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p)) - eold
        dg = deltag_growth(solubd, nso+dn, nso)
        lpref = lnprefsol(-dn, solution%monos, solution%vol)
        !lpref = lpref + log(totally_not_kinetic_factor(solubd, number_monos(p)))

        ret = acceptsornot(de, lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = solution%monos - dn
        else
          call sac_change_part_size(st, ibox, p, nso)
        end if
        ! only if acc ?
    
        deltaForm_precip(f) = deltaForm_precip(f) + (dg-lpref)/dn
        deltaE_precip(f) = deltaE_precip(f) + de%ep/dn
        deltaU_precip_count(f) = deltaU_precip_count(f) + 1
      else
      !   ! new particle
        p = st%boxes(ibox)%n+1
        addordelpart = .true.
        do k=1,3
          pos(k) = ran2(st%rng)
        end do
        part = vsgc_make_part(st, pos, nmonos_to_rad(solubd, dn))
        ! call test_nrg(st, ibox, "Before Precip Add")

        ! print *, eold%ep, " ", st%boxes(ibox)%ener%ep, st%boxes(ibox)%n
        call sac_add_part_box(st, ibox, part, dn)
        de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        ! print *, "aft: ", eold%ep+de%ep

        dg = deltag_part(solubd, dn, lbox%volume)
        lpref = lnprefsol(-dn, solution%monos, solution%vol)

        ret = acceptsornot(de, lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = solution%monos - dn
          part_created_count = part_created_count + 1_8
        else
          call sac_del_part_box(st, ibox, p)
        end if
      end if
    else
      ret = mv_cir_rej
    end if
  end if


  end associate
contains

 

  integer function acceptsornot(de, dp) result(rep)
    type(NRGRec), intent(inout) :: de
    real(8), intent(in) :: dp
    if(de%noverlap) then
      if(metropolis(de%ep/kT - dp, st%rng)) then
        rep = mv_acc
        st%boxes(ibox)%ener%ep = st%boxes(ibox)%ener%ep + de%ep
      else
        rep = mv_nrj_rej
      end if
    else
      rep = mv_hs_rej
    end if
  end function
    
end function

integer function particle_mono_exchange_vkin(st, nr, dissolution) result(ret)
  use iso_fortran_env
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  logical, intent(in) :: dissolution
  integer :: ibox, p,k, nsmall, j, i, np_tot
  integer(8) :: dn, m0, m1, nso
  real(8) :: r1p, r2p, fac, r1pold, r2pold, q1old, q2old, q1p, &
    q2p, pos(3), vpo, vpn, lpref, dg
  type(NRGRec) :: eold, eold2, de
  real(8), parameter :: onethrd = 1._8/3._8
  type(Particule) :: part
  logical :: addordelpart
  integer, save :: counter = 0
  counter = counter + 1
  ret= mv_acc
  ibox = 1 ! ⚠ for now
  associate(lbox => st%boxes(ibox))
  addordelpart = .false.
  np_tot = st%boxes(ibox)%n
  
  ! selects a particule
  p = floor(ran2(st%rng)*np_tot) + 1

  solubd%kinetic_constant = 20./mean_particle_radius(st, ibox, r2p)**2
  dn = 1_8 + floor(ran2(st%rng)*totally_not_kinetic_factor(solubd, number_monos(p)))
  if(dissolution) then
    ! Dissolution
    if(np_tot > 0) then
      if(number_monos(p) >= dn + minimum_monos_particles) then
        dg = deltag_growth(solubd, number_monos(p)-dn, number_monos(p))
        lpref = lnprefsol(dn, solution%monos, solution%vol)
        if(bulk_exchange_create_destroy) lpref = lpref + log(real(lbox%n,8)/real(lbox%n+1, 8))

        nso = number_monos(p)
        eold = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        call sac_change_part_size(st, ibox, p, nso-dn)
        de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p)) - eold

        ret = acceptsornot(de, lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = dn + solution%monos
        else
          call sac_change_part_size(st, ibox, p, nso)
        end if
      else
        ret = mv_cir_rej
      end if
    else
      ret = mv_cir_rej
    end if
  else
    ! Growth
    if(solution%monos >= dn) then
      if(p <= np_tot) then

        nso = number_monos(p)
        eold = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        call sac_change_part_size(st, ibox, p, nso+dn)
        de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p)) - eold
        dg = deltag_growth(solubd, nso+dn, nso)
        lpref = lnprefsol(-dn, solution%monos, solution%vol)

        ret = acceptsornot(de, lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = solution%monos - dn
        else
          call sac_change_part_size(st, ibox, p, nso)
        end if
      end if
    else
      ret = mv_cir_rej
    end if
  end if

  end associate
contains

  integer function acceptsornot(de, dp) result(rep)
    type(NRGRec), intent(inout) :: de
    real(8), intent(in) :: dp
    if(de%noverlap) then
      if(metropolis(de%ep/kT - dp, st%rng)) then
        rep = mv_acc
        st%boxes(ibox)%ener%ep = st%boxes(ibox)%ener%ep + de%ep
      else
        rep = mv_nrj_rej
      end if
    else
      rep = mv_hs_rej
    end if
  end function
    
end function


integer function particle_bulk_dissolution(st, nr) result(ret)
  use iso_fortran_env
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  ret = particle_bulk_reformation(st, nr, .true.)
end function

integer function particle_bulk_precipitation(st, nr) result(ret)
  use iso_fortran_env
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  ret = particle_bulk_reformation(st, nr, .false.)
end function

integer function particle_bulk_reformation(st, nr, dissolution) result(ret)
  use iso_fortran_env
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  logical, intent(in) :: dissolution
  integer :: ibox, p,k, nsmall, j, i, ismall, psmall
  integer(8) :: dn, m0, m1, d_dn
  real(8) :: r1p, r2p, fac, r1pold, r2pold, q1old, q2old, q1p, &
    q2p, pos(3), vpo, vpn, lpref, dg
  type(NRGRec) :: eold1, eold2, de
  real(8), parameter :: onethrd = 1._8/3._8
  type(Particule) :: part
  logical :: addordelpart
  ret= mv_acc
  
  ibox = 1 ! ⚠ for now
  associate(lbox => st%boxes(ibox))
  m0 = solution%monos + sum(number_monos(1:lbox%n))

  nsmall = 0
  do i=1,st%boxes(ibox)%n
    if(number_monos(i) <= maximum_monos_particles) nsmall = nsmall + 1
  end do
  d_dn = maximum_monos_particles - minimum_monos_particles
  if(dissolution) then
    ! Dissolution
    ! selects a particule whose nmono <= maximum_monos_particles
    if(nsmall > 0) then
      psmall = floor(ran2(st%rng)*nsmall) + 1
      ismall = 0
      do p=1, st%boxes(ibox)%n
        if(number_monos(p) <= maximum_monos_particles) then
          ismall = ismall + 1
          if(ismall == psmall) exit
        end if
      end do
  
      dn = number_monos(p)
      if(dn<=maximum_monos_particles) then
        ! delete the particle
        dg = - deltag_part(solubd, dn, lbox%volume)
        lpref = lnprefsol(dn, solution%monos, solution%vol)
        part = lbox%parts(p)
        de = - nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        call sac_del_part_box(st, ibox, p)
        ret = acceptsornot(de, log(real(nsmall,8)/d_dn) + lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = dn + solution%monos
          part_deleted_count_integ = part_deleted_count_integ + 1
        else
          call sac_add_part_box(st, ibox, part, dn)
        end if
      else
      ret = mv_cir_rej    
      end if
    else
      ret = mv_cir_rej
    end if
  else
    ! Precipitation
    dn = minimum_monos_particles+ floor(ran2(st%rng) * (d_dn))
    if(solution%monos >= dn) then
      p = st%boxes(ibox)%n+1
      if(.true.) then
        ! new particle
        do k=1,3
          pos(k) = ran2(st%rng)
        end do
        part = vsgc_make_part(st, pos, nmonos_to_rad(solubd, dn))
        call sac_add_part_box(st, ibox, part, dn)
        !call change_part_radius(st, ibox, p, nmonos_to_rad(solubd, dn))
        dg = deltag_part(solubd, dn, lbox%volume)
        lpref = lnprefsol(-dn, solution%monos, solution%vol)
        de = nr%energy_1p(st,st%boxes(ibox),st%boxes(ibox)%parts(p))
        ret = acceptsornot(de, log(real(d_dn,8)/(nsmall+1)) + lpref-dg)
        if(ret == mv_acc) then
          if(finite_monomer_reservoir) solution%monos = solution%monos - dn
          part_created_count_integ = part_created_count_integ + 1_8
        else
          call sac_del_part_box(st, ibox, p)
        end if
      ! else
      !   ret = mv_cir_rej
      end if
    else
      ret = mv_cir_rej
    end if
  end if
  end associate
contains

  integer function acceptsornot(de, dp) result(rep)
    type(NRGRec), intent(inout) :: de
    real(8), intent(in) :: dp
    if(de%noverlap) then
      if(metropolis(de%ep/kT - dp, st%rng)) then
        rep = mv_acc
        st%boxes(ibox)%ener%ep = st%boxes(ibox)%ener%ep + de%ep
      else
        rep = mv_nrj_rej
      end if
    else
      rep = mv_hs_rej
    end if
  end function
    
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
    if(npt_count_disappeared) then
      truen = labox%n - disappeared_particles(st, ibox)
    else
      truen = labox%n
    end if

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

  integer :: i, fo
  real(8) :: dr

  call histor8_set_distribution(p)

  dr = 1._8 / p%dri
  do i=lbound(p%h,1),ubound(p%h,1)
    write(unit, '(2E16.7)') p%rmin + (i-0.5_8)*dr, p%h(i)
  end do

end subroutine


end module

module inputp
!! Definition of input parameters.
implicit none

  type Input
    !! Type containing input parameters.
    character(len=:), allocatable :: file_prefix
    character(len=:), allocatable :: fich_table
    character(len=:), allocatable :: asimu
    character(len=:), allocatable :: smod_pot
    character(len=:), allocatable :: smod_box
    character(len=:), allocatable :: smod_alg
    logical :: old_fashion = .false.
    integer :: step_w = -1
    integer :: step_nrj = -1
    integer :: per_reuss_trans = -1
    integer :: seed = 666
    integer :: max_number_of_particles = -1
    integer :: number_of_boxes = 1
    real(8) :: dr_press = 0.05
    real(8) :: temp = 300.d0
    real(8) :: debye_length = 4.32d0
    real(8) :: length_prefactor = 0.715d0
    real(8) :: eps = 0.001d0
    real(8) :: reuss_trans_ref = 0.5
    real(8) :: trans_ampl = 3.
    real(8) :: rcutoff = -1.d0
    real(8) :: trans_ampl_wall = -1.d0
    real(8) :: boxlx = 1.d0
    real(8) :: boxly = 1.d0
    real(8) :: boxlz = 1.d0
    real(8) :: boxalpha = 90.d0
    real(8) :: boxbeta = 90.d0
    real(8) :: boxgamma = 90.d0
    real(8) :: phi_inp=-1.D0
    real(8) :: delta_nrj_tol = 1.D-5
    real(8) :: potential_rmax = 1.d99
    integer :: initial_step = 0
    real(8) :: dr_calc_cutoff = 0.1
    logical :: prod = .true.
    logical :: conf_temp_on = .false.
    logical :: std_gibbs_volex = .false.
    integer :: maxcycle_stab = 0
    integer :: maxcycle_calc = 1
    integer :: step_job = 0
    integer :: step_log = 0
    real(8) :: pref
    integer :: per_reuss_vol
    integer :: stepfluctu = -1
    integer :: step_bigs = -1
    integer :: ntrials_rosbth = 10
    integer :: n_rospechmax = -1
    real(8) :: dr_rdf = 0.1
    real(8) :: dq = 0.002
    real(8) :: prob_swap = 0.d0
    real(8) :: reuss_vol_ref
    real(8) :: prob_swap_box = 0.d0
    real(8) :: prob_box2box = 0.d0
    real(8) :: prob_volex = 0.d0
    real(8) :: dv_volex = 0.d0
    real(8) :: prob_volparex = 0.d0
    real(8) :: prob_b2b_rosbth = 0.d0
    real(8) :: prob_swap_impur = 0.d0
    real(8) :: prob_swap_impur_rosbth = 0.d0
    real(8) :: prob_swap_impur_rosbth2 = 0.d0
    real(8) :: prob_volmod = 0.d0
    real(8) :: dv_vol = 0.d0
    real(8) :: rp_impur_rosbth=0.5
    real(8) :: rvois_impur_rosbth=-1.D0
    logical :: statcell=.false.
    logical :: chempwid_doit = .false.
    integer :: pZdens = -1
    logical :: fluctudens=.false.
    logical :: nrg_ergo_analysis=.false.
    logical :: use_bin_traj=.false.
    character(len=:), allocatable :: fich_fluctu
    integer :: widom_period = -1
    integer :: widom_samples = 1
    integer :: seed_analysis = 16385673
    integer :: internal_steps = -1
    logical :: semigrand = .false.
  end type Input


contains

  subroutine init_input(inp)
    type(Input) :: inp
    inp%file_prefix = "pcmc"  
    inp%fich_table = "table.dat"
    inp%asimu = "NVT"
    inp%smod_pot = "yukawa"
    inp%smod_box = "generic"
    inp%smod_alg = "standard"
    inp%fich_fluctu = "fluctuations_input.txt"
  end subroutine

  logical function input_find_keyval(inp,key, val)
    !! Parse (*key*-*val*) pairs to corresponding field.
    !! Returns `.true.` if key is found
    use readparam, only: read_logical
    type(Input) :: inp
    character(*), intent(in) :: key, val
    input_find_keyval = .true.

    select case(key)
    case("Bjerrum_length")
      read(val,*) inp%length_prefactor
    case("temperature") 
      read(val,*) inp%temp
    case("Debye_length") 
      read(val,*) inp%debye_length
    case("random_seed") 
      read(val,*) inp%seed
    case("simulation_type") 
      inp%asimu = val
    case("snapshot_period") 
      read(val,*) inp%step_w
    case("energy_output_period") 
      read(val,*) inp%step_nrj
    case("translation_update_period") 
      read(val,*) inp%per_reuss_trans
    case("translation_interval")
      read(val,*) inp%trans_ampl
    case("translation_acceptance") 
      read(val,*) inp%reuss_trans_ref
    case("hsp_interval") 
      read(val,*) inp%dr_press
    case("cutoff_estimation_tolerance") 
      read(val,*) inp%eps
    case("table_file") 
      inp%fich_table = val
    case("old_charges") 
      inp%old_fashion = read_logical(val)
    case("cutoff_radius") 
      read(val,*) inp%rcutoff
    case("box_length_X") 
      read(val,*) inp%boxLx
    case("box_length_Y") 
      read(val,*) inp%boxLy
    case("box_length_Z") 
      read(val,*) inp%boxLz
    case("alpha_box_angle") 
      read(val,*) inp%boxalpha
    case("beta_box_angle") 
      read(val,*) inp%boxbeta
    case("gamma_box_angle") 
      read(val,*) inp%boxgamma
    case("number_of_boxes") 
      read(val,*) inp%number_of_boxes
    case("max_number_of_particles") 
      read(val,*) inp%max_number_of_particles
    case("box_type") 
      inp%smod_box = val
    case("potential_type") 
      inp%smod_pot = val
    case("energy_algorithm") 
      inp%smod_alg = val
    case("volume_fraction") 
      read(val,*) inp%phi_inp
    case("delta_energy_tolerance") 
      read(val,*) inp%delta_nrj_tol
    case("initial_step") 
      read(val,*) inp%initial_step
    case("translation_interval_wall") 
      read(val,*) inp%trans_ampl_wall
    case("output_prefix") 
      inp%file_prefix = val
    case("cutoff_estimation_precision") 
      read(val,*) inp%dr_calc_cutoff
    case("production") 
      inp%prod  = read_logical(val)
    case("equilibration_steps") 
      read(val,*) inp%maxcycle_stab
    case("number_of_steps")
      read(val,*) inp%maxcycle_calc
    case("calculation_period") 
      read(val,*) inp%step_job
    case("rdf_bin_length") 
      read(val,*) inp%dr_rdf
    case("wave_number_interval") 
      read(val,*) inp%dq
    case("swap_probability") 
      read(val,*) inp%prob_swap
    case("reference_pressure") 
      read(val,*) inp%Pref  
    case("fluctuations_period") 
      read(val,*) inp%stepfluctu
    case("configurational_temperature") 
      inp%conf_temp_on = read_logical(val)
    case("log_period") 
      read(val,*) inp%step_log
    case("full_structure_period") 
      read(val,*) inp%step_bigs
    case("interbox_swap_probability") 
      read(val,*) inp%prob_swap_box
    case("box_to_box_probability") 
      read(val,*) inp%prob_box2box
    case("volume_exchange_probability") 
      read(val,*) inp%prob_volex
    case("exchange_deltavolume") 
      read(val,*) inp%dv_volex
    case("volume_particle_exchange_probability") 
      read(val,*) inp%prob_volparex
    case("standard_volume_change") 
      inp%std_gibbs_volex = read_logical(val)
    case("on_the_fly_chemical_potential") 
      inp%chempwid_doit = read_logical(val)
    case("particle_exchange_rosenbluth_probability") 
      read(val,*) inp%prob_b2b_rosbth
    case("number_of_rosenbluth_trials") 
      read(val,*) inp%ntrials_rosbth
    case("maximum_family_number_rosenbluth") 
      read(val,*) inp%n_rospechmax
    case("impure_interbox_swap_probability") 
      read(val,*) inp%prob_swap_impur
    case("cell_statistics") 
      inp%statcell = read_logical(val)
    case("interbox_swap_impure_rosbth_probability") 
      read(val,*) inp%prob_swap_impur_rosbth
    case("interbox_swap_impure_rosbth_radius") 
      read(val,*) inp%rp_impur_rosbth
    case("interbox_swap_impure_rosbth_2_probability")  
      read(val,*) inp%prob_swap_impur_rosbth2
    case("z_density_period")
      read(val,*) inp%pZdens
    case("volume_change_probability") 
      read(val,*) inp%prob_volmod
    case("deltavolume")
      read(val,*) inp%dv_vol
    case("density_fluctuations_file") 
      read(val,*) inp%fich_fluctu
    case("density_fluctuations") 
      inp%fluctudens = read_logical(val)
    case("potential_max_separation")
      read(val,*) inp%potential_rmax
    case("energy_ergo_analysis")
      inp%nrg_ergo_analysis = read_logical(val)
    case("use_binary_trajectory")
      inp%use_bin_traj = read_logical(val)
    case("widom_period")
      read(val, *) inp%widom_period
    case("widom_samples")
      read(val, *) inp%widom_samples
    case("internal_steps")
      read(val, *) inp%internal_steps
    case("semigrand")
      inp%semigrand = read_logical(val)
    case default
      input_find_keyval = .false.        
    end select
  end function input_find_keyval

  subroutine config_to_input(list, inp)
    !! Sets `inp` from key-val pair list `list`
    use readparam
    type(Input) :: inp
    type(ParamList) :: list
    character(len=:), allocatable :: key,val
    integer :: keyLen
    call start_iter_keyval(list)
    do while(iter_keyval_in(list))
      call get_current_keyval(list, key, val)
      keyLen = len(key)
      if(input_find_keyval(inp,key,val)) then
        call access_current(list)
      end if
      call param_list_next(list)
    end do
  end subroutine config_to_input

end module inputp
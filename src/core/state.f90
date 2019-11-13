module state
  !! Defines the base types and the state type.
  use random2
  use geom
  use readparam, only: ParamList
  use celldec
  use inputp
  implicit none

  type Distrib
    !! Distribution
    integer :: nfam  !! number of families
    integer,allocatable :: pop(:)    !! distribution populations array
    real(8),allocatable :: rad(:)    !! distribution radii array
    real(8),allocatable :: charge(:) !! distribution charges array
  end type Distrib
  
  type NrgRec
    !! Energy
    real(8) :: Ep !! Potential energy
    logical :: noverlap !! false if there is no overlap
  end type NrgRec

  type Pressure
    !! Pressure
    real(8) :: virp !! Virial pressure
    real(8) :: hsp !! Hard sphere pressure
  end type Pressure

  interface operator(+)
    module procedure NrgRec_add
  end interface
  interface operator(-)
    module procedure NrgRec_sub
  end interface
  interface operator(-)
    module procedure NrgRec_neg
  end interface

  type Box
     !! Box structure
     integer :: num   !! number (id)
     integer :: n    !! number of particles inside
     real(8) :: tens(3,3) !! box tensor 
     real(8) :: met(3,3)  !! box metric tensor
     real(8) :: volume    !! volume
     integer, allocatable :: fam(:)  !! array of family populations
     type(Particule), allocatable :: parts(:)  !! array of particles
     type(NrgRec) :: Ener  !! Energy
     type(Pressure) :: press  !! pressures
  end type Box


  type BoxRes
    !! Average quantities of a box
    real(8) :: volm  !! Mean volume
    real(8) :: vol2m !! volume  variance
    real(8) :: nm    !! mean number of particles
    real(8) :: ep    !! mean potential energy
    real(8) :: p     !! mean pressure
    real(8) :: virp  !! mean virial pressure
    real(8) :: hsp   !! mean hard sphere pressure
    real(8) :: var_ep   !! potential energy variance
    real(8) :: var_hsp !! hard spere pressure variance
    real(8) :: var_virp !! virial pressure variance
    real(8) :: e2c  !! total energy variance
    real(8) :: e    !! total energy
    integer :: nbetamu   !! ?
    real(8),allocatable :: popm(:) !! mean population per family
    real(8) :: hetv     !!  volume fraction
  end type BoxRes

  type Pmconf
    !! Trajectory file handle
    integer :: unit
    integer :: kind
  end type Pmconf
  
  type StepData
    !! Counters updated at each MC step.
    integer :: nnn
    integer(8) :: reuss_trans,trans_eff,reuss_trans_p,trans_eff_p,reuss_trans_outwall_p,reuss_trans_wall_p,trans_outwall_p,&
        trans_wall_p
    integer(8) :: step_count,step_nrj_count, step_calc_n
    integer :: step_nrj_n,step_w_n,step_uptr
    logical :: cycle_update_trans
    integer :: cycle_first_step
  end type StepData

  type RtParams
    real(8) :: dx_trans
    real(8) :: dx_trans_wall
    real(8) :: Lmin
  end type RtParams

  type NrgData
    !! Used by energy calculations
    integer :: mode_combo,mode_pot,mode_box,mode_dec
    real(8) :: inv_debye_length, elpot_pref
    real(8) :: kappa2, rcutoff_max
    real(8) :: inv_rc2, minus2_inv_rc
    real(8), allocatable :: eff_charge(:)
  end type NrgData

  integer, parameter :: max_boxes = 32 !! Maximum number of boxes. Should be Enough™
  type PcmcState
    !! Represents the current state of a simulation.
    type(Input) :: inp                   !! input parameter
    integer :: ntotbox
    type(Box) :: boxes(max_boxes)        !! Array of all boxes
    type(BoxRes) :: res(max_boxes)      !! Array of intermediate results for each box
    type(Distrib) :: dist                !! Distribution
    integer :: step
    type(Pmconf) :: finst
    type(StepData) :: sdata
    type(Ran2State) :: rng
    type(Ran2State) :: rng_an   ! for analysis
    type(RtParams) :: rpar
    type(NrgData) :: ndat
    type(ParamList) :: config
    ! some other vars
    integer :: Nbuff=-1
    real(8) :: rcrc,rc,lmin,lmax
    real(8) :: rcut_ap = 0.d0
    integer :: transss, reuss_trans_outwall, reuss_trans_wall
    integer :: trans_wall=0,trans_outwall=0
    ! Modes
    integer :: mode_box,mode_pot,mode_dec,mode_combo
    logical :: celllist
    type(CellDecomp) :: clistdecomp(max_boxes)
  end type PcmcState
  
  real(8),parameter :: augment_dx=1.052631579D0          !! factor to increase translation amplitude
  real(8),parameter :: reduit_dx=0.95D0                  !! factor to decrease translation amplitude
  real(8), parameter :: dv_volex_max = 10.d0             !! maximum ΔV for volume changes.

  real(8),parameter :: PI=3.1415926535897931D0
  real(8),parameter :: DEG2RAD=1.74532925199432955D-002  !! Degress to radiant conversion factor
  real(8),parameter :: k_B = 1.38064852d-23              !! Boltzmann constant 
  real(8),parameter :: N_A = 6.022140857d23              !! Avogadro number

  ! Units
  real(8), parameter :: kT=1.d0                          !! value of k_B T. Set to 1. because energy is in k_B T units.
  real(8), parameter :: unitl=1.0d-9                     !! Length unit in m. 
  real(8), parameter :: mBarSi=100.d0                    !! 1 millibar in Pa.
  real(8), parameter :: kJpMolpK=N_A*k_B/1000            !! kJ/mol/K
  real(8), parameter :: unitEsT=kJpMolpK                 !! Energy convertion factor.
  real(8), parameter :: unitPsT=k_B/unitL**3             !! Pressure convertion factor from (T units)/nm³ to Pa. 


  ! For trajectories
  integer, private :: bconfMag = int(z'0000BC0F')
  integer, parameter, private :: bconfHead = 100007, bconfBox = 1000008, bconfBoxEnd = 1000080, &
    bconfFam = 1000009, bconfPos = 10000010, bconfRad = 10000011, huge4 = huge(4)

  integer,parameter :: ckConf=1  !!  ascii inst kind
  integer,parameter :: ckBConf=2  !!  binary inst kind

  private :: write_conf_frame, read_conf_frame, write_bconf_frame, read_bconf_frame

contains

  subroutine box_set_tensor(b, p)
    type(Box) :: b
    real(8), intent(in) :: p(3,3)
    b%tens = p
    call calc_metric(p, b%met)
    b%volume = volume_box(P)
  end subroutine box_set_tensor

  type(NrgRec) function NrgRec_add(e1,e2) result(res)
    !! Adds two energies. Overlap state is propagated.
    type(NrgRec), intent(in) :: e1, e2
    res%noverlap = e1%noverlap .and. e2%noverlap
    res%ep = e1%ep + e2%ep
  end function NrgRec_add

  type(NrgRec) function NrgRec_sub(e1,e2) result(res)
    !! Substracts two energies. Overlap state is propagated.
    type(NrgRec), intent(in) :: e1, e2
    res%noverlap = e1%noverlap .and. e2%noverlap
    res%ep = e1%ep - e2%ep
  end function NrgRec_sub

  type(NrgRec) function NrgRec_neg(e1) result(res)
    !! Negates one energy. Overlap state is propagated.
    type(NrgRec), intent(in) :: e1
    res%noverlap = e1%noverlap
    res%ep = - e1%ep
  end function NrgRec_neg

  subroutine deallocate_boxes(st)
    type(PcmcState) :: st
    integer :: k
    do k=1,st%ntotbox
      st%boxes(k)%num = 0
      st%boxes(k)%n = 0
      deallocate(st%boxes(k)%parts)
      if(st%dist%nfam > 0) deallocate(st%boxes(k)%fam)
    end do
  end subroutine deallocate_boxes

  subroutine allocate_boxes(st,n_buff)
    type(PcmcState) :: st
    integer,intent(in) :: n_buff
    integer :: k
    do k=1,st%ntotbox
      st%boxes(k)%num = k
      allocate(st%boxes(k)%parts(n_buff))
      if(st%dist%nfam > 0) allocate(st%boxes(k)%fam(st%dist%nfam))
    end do
  end subroutine allocate_boxes

  subroutine allocate_multires(st)
    integer :: k
    type(PcmcState) :: st
    do k=1,st%ntotbox            
      allocate(st%res(k)%popm(st%dist%nfam))
    end do
  end subroutine allocate_multires

  subroutine deallocate_multires(st)
    type(PcmcState) :: st
    integer :: k
    do k=1,st%ntotbox
      deallocate(st%res(k)%popm)
    end do
  end subroutine deallocate_multires

  logical function read_distrib(dist,distfile)
    implicit none
    type(Distrib) :: dist
    character(*),intent(in) :: distfile
    logical :: estla
    integer :: i, uf
    inquire(file=distfile,exist=estla)
    if(.not.estla) then
      read_distrib = .false.
      return
    endif
    open(newunit=uf,file=distfile,status='old')
    read(uf,*) dist%nfam
    allocate(dist%pop(dist%nfam),dist%rad(dist%nfam),dist%charge(dist%nfam))
    do i=1,dist%nfam
      read(uf,*) dist%pop(i),dist%rad(i),dist%charge(i)
    enddo
    close(uf)
    read_distrib=.true.
  end function read_distrib

  subroutine make_distrib_from_conf(st)
    type(PcmcState) :: st
    integer :: i,k, nr, fam
    nr=1
    do k=1,st%ntotbox
      do i=1,st%boxes(k)%n
        nr = max(st%boxes(k)%parts(i)%famille, nr)
      enddo
    enddo

    if( allocated(st%dist%pop)) then
      deallocate(st%dist%pop, st%dist%rad, st%dist%charge)
    endif
    st%dist%nfam = max(nr, st%dist%nfam)
    allocate(st%dist%pop(st%dist%nfam), st%dist%rad(st%dist%nfam), st%dist%charge(st%dist%nfam))
    st%dist%pop = 0
    do k=1,st%ntotbox
      do i=1,st%boxes(k)%n
        fam = st%boxes(k)%parts(i)%famille
        st%dist%pop(fam) = st%dist%pop(fam) + 1
        st%dist%rad(fam) = st%boxes(k)%parts(i)%rayon
      enddo
    enddo
  end subroutine make_distrib_from_conf
  
  subroutine part_set_props(part, st)
    type(PcmcState), intent(in) :: st
    type(Particule) :: part
    integer :: fam
    fam = part%famille
    !part%rayon = st%dist%rad(fam)
    part%ech = st%ndat%eff_charge(fam)
  end subroutine part_set_props
  
  subroutine affect_radius_from_distrib(st)
    type(PcmcState) :: st
    integer :: k,i
    real(8) :: kp
    kp = 1.d0 / st%inp%debye_length
    do k=1,st%ntotbox
      do i=1, st%boxes(k)%n
        st%boxes(k)%parts(i)%rayon = st%dist%rad(st%boxes(k)%parts(i)%famille)
      end do
    end do
  end subroutine affect_radius_from_distrib

  subroutine update_box_fam(b,dist)
    implicit none
    type(Box) :: b
    type(Distrib), intent(in) :: dist
    integer :: i,j
    do j=1,dist%nfam
      b%fam(j) = 0
    end do
    do i=1,b%n
      j=b%parts(i)%famille
      b%fam(j) = b%fam(j) + 1
    end do
  end subroutine update_box_fam

  integer function get_mode_box(s_mod) result(m)
    implicit none
    character(*),intent(in) :: s_mod

    select case(s_mod)
    case("generic")
      m = 0
    case("cubic")
      m = 1
    case("cubicXY")
      m = 2
    case("hexagonal")
      m = 3
    case default
      m = 0
    end select
  end function get_mode_box

  integer function get_mode_alg(s_mod) result(m)
    implicit none
    character(*),intent(in) :: s_mod

    select case(s_mod)
    case("standard")
      m = 0
    case("celldec")
      m = 1
    case default
      m = 0
    end select
  end function get_mode_alg

  ! Initialize box and nrg alg. modes
  subroutine init_modes(st)
    type(PcmcState) :: st

    st%mode_pot = 0 ! Actual value set at nrg_init()
    st%mode_box = get_mode_box(st%inp%smod_box)
    st%mode_dec = get_mode_alg(st%inp%smod_alg)

    if(st%mode_dec == 1) then
      st%celllist = .true.
    end if

  end subroutine init_modes

  subroutine init_distribution(dist,fich,status)
    type(Distrib) :: dist
    character(*) :: fich
    integer :: status
    logical :: inestla
    integer :: i, uf
    status = 0
    inquire(file=fich,exist=inestla)
    if(.not.inestla) then
      status = -1
      return
    endif
    open(newunit=uf,file=fich,status='old')
    read(uf,*) dist%nfam
    allocate(dist%pop(dist%nfam),dist%rad(dist%nfam),dist%charge(dist%nfam))
    do i=1,dist%nfam
      read(uf,*) dist%pop(i),dist%rad(i),dist%charge(i)
    enddo
    close(uf)
  end subroutine init_distribution

  subroutine read_cfg_file(st,fichpar)
    use readparam
    implicit none
    type(PcmcState) :: st
    character(*), intent(in) :: fichpar
    integer :: cfu
    open(newunit=cfu, file=fichpar, status='old')
    call init_param_list(st%config)
    if(read_configuration(cfu, st%config) /= 0) then
      write(error_unit, '(A)') "Error in config file '", trim(fichpar), "'."
      stop 1
    end if
    call config_to_input(st%config, st%inp)
    call configuration_finish(st)
    close(cfu)
  end subroutine read_cfg_file

  subroutine configuration_finish(st)
    type(PcmcState) :: st
    if(st%inp%number_of_boxes < 1) then
      st%ntotbox = st%inp%number_of_boxes
    else if(st%inp%number_of_boxes > max_boxes) then
      st%ntotbox = max_boxes
    else
      st%ntotbox = st%inp%number_of_boxes
    end if
    if(st%Nbuff < st%inp%max_number_of_particles) then
      st%Nbuff = st%inp%max_number_of_particles
    endif
    st%rpar%dx_trans = st%inp%trans_ampl
    st%rpar%dx_trans_wall = st%inp%trans_ampl_wall
  end subroutine configuration_finish

  subroutine init_box_from_configuration(st,boxLx,boxLy,boxLz,boxalpha,boxbeta,boxgamma,phi, unit)
    use iso_fortran_env
    implicit none
    type(PcmcState) :: st
    real(8) ::  boxLx,boxLy,boxLz,boxalpha,boxbeta,boxgamma
    real(8), intent(in), optional :: phi
    integer, intent(in), optional :: unit
    real(8) :: volp, dump1,volume
    real(8) :: boxtens(3,3)
    integer :: i, logu
    if(present(unit)) then
      logu = unit
    else
      logu = output_unit
    end if
    if(present(phi) .and. phi > 0.d0) then
      call create_box_from_3l3a(boxtens,boxLx,boxLy,boxLz,boxalpha*DEG2RAD,boxbeta*DEG2RAD,boxgamma*DEG2RAD)
      volume = volume_box(boxtens)
      volp = 4.d0/3.d0*pi*sum(st%dist%pop*st%dist%rad**3)
      dump1= volp/(st%ntotbox*volume)
      write(logu,'(A,G11.4)') "Volume  fraction from box parameters: ", volp/volume
      write(logu,'(A,G11.4)') "Setting to requested volume fraction: ", phi
      boxLx = boxLx*(dump1/phi)**(1./3.)
      boxLy = boxLy*(dump1/phi)**(1./3.)
      boxLz = boxLz*(dump1/phi)**(1./3.)
    end if
    if(st%ntotbox<1) then
      write(error_unit,*) "Invalid number of box : ",st%ntotbox,"."
      stop 6
    end if
    call allocate_boxes(st,st%Nbuff)
    do i=1,st%ntotbox
      call create_box_from_3l3a(st%boxes(i)%tens,boxLx,boxLy,boxLz,boxalpha*DEG2RAD,boxbeta*DEG2RAD,boxgamma*DEG2RAD)
      call calc_metric(st%boxes(i)%tens, st%boxes(i)%met)
    enddo
  end subroutine init_box_from_configuration

subroutine close_pmconf(f)
  type(Pmconf) :: f
  close(f%unit)
end subroutine

function open_pmconf(path, old, bin)
  type(Pmconf) :: open_pmconf
  integer :: fu
  character(*) :: path
  logical, optional :: old, bin
  logical :: tbin, told
  integer :: hid  
  if(present(bin)) then
    tbin = bin
  else
    tbin = .false.
  end if
  if(present(old)) then
    told = old
  else
    told = .false.
  end if

  ! Guess the format
  if(told) then
    open(newunit=fu, file=path, form='unformatted', access='stream', status='old', action='read')
    read(fu) hid
    if(hid == bconfMag) then
      tbin = .true.
    else
      tbin = .false.
    end if
    close(fu)
  end if
  if(tbin) then
    open_pmconf%kind = ckBConf
    if(told) then
      open(newunit=fu, file=path, status='old', form='unformatted', access='stream',action='read')
    else
      open(newunit=fu, file=path, status='replace', form='unformatted', access='stream', action='write')
    end if
  else
    open_pmconf%kind = ckConf
    if(told) then
      open(newunit=fu, file=path, status='old', action='read')
    else
      open(newunit=fu, file=path, status='replace', action='write')
    end if
  end if
  open_pmconf%unit = fu
end function

  !! Write the current configuration to the file unit fu.
  !! If optional write_cfg is set to .false., do not write some parameters.
  subroutine write_pmconf_frame(st,f,no_cfg)
    implicit none
    type(PcmcState) :: st
    type(Pmconf) :: f
    logical, optional :: no_cfg
    logical :: nocfg
   
    if(present(no_cfg)) then
      nocfg = no_cfg
    else
      nocfg = .false.
    end if

    select case(f%kind)
    case(ckConf)
      call write_conf_frame(st,f%unit,nocfg)
      flush f%unit
    case(ckBConf)
      call write_bconf_frame(st,f%unit)
      flush f%unit
    case default
      stop 666
    end select


   end subroutine write_pmconf_frame

  !! Read the current configuration from the file unit fu.
  !! If realloc is set to `.true.` the boxes are reallocated.
  !! iostat for I-O error handling.
  subroutine read_pmconf_frame(st,f,realloc, iostat)
    implicit none
    type(PcmcState) :: st
    type(Pmconf) :: f
    integer,intent(out)::  iostat
    logical,intent(in) :: realloc
    iostat = -10 !
    if( allocated(st%boxes(1)%parts)) then
      if(realloc) then
        call deallocate_boxes(st)
      end if
    end if
    if(.not. allocated(st%boxes(1)%parts)) then
      call allocate_boxes(st,st%Nbuff)
    end if

    select case(f%kind)
    case(ckConf)
      call read_conf_frame(st,f%unit, iostat)
    case(ckBConf)
      call read_bconf_frame(st,f%unit, iostat)
    case default
      stop 666
    end select

  end subroutine read_pmconf_frame

  ! text variant
  subroutine write_conf_frame(st,fu,no_cfg)
    use geom
    implicit none
    type(PcmcState) :: st
    integer,intent(in) :: fu
    logical, optional :: no_cfg
    logical :: wcfg
   
    integer :: i,j,k
    character(6) :: name
    real(8) :: vec3(3)
    name="      "
    if(present(no_cfg)) then
      wcfg = .not. no_cfg
    else
      wcfg = .true.
    end if

    write(fu,'("step ",I8)') st%step
    write(fu,'("ntotbox =",I8)') st%ntotbox
    if(wcfg) then
      write(fu, '("#$config")')
      call write_changed_paras(st,fu)
      write(fu, '("#$endconfig")')
    else
      write(fu,*)
    end if
    do k=1,st%ntotbox
      write(fu,'("boxtensor = ")')
      do i=1,3
        write(fu,'(3ES16.7)') st%boxes(k)%tens(i,1:3)
        
      end do
      j=st%boxes(k)%n
      write(fu,'("number of particles =",I8)') j
      do i=1,j
        vec3=matvec3_prod(st%boxes(k)%tens,st%boxes(k)%Parts(i)%pos)
        !write(name,'(a3,i0)') "SPH",st%boxes(k)%Parts(i)%famille
        write(fu,'(3ES16.7,I6,ES16.7)') vec3(1:3),&
              st%boxes(k)%Parts(i)%famille,st%boxes(k)%Parts(i)%rayon
      end do
      write(fu,*)
    end do
    write(fu,'(A)') "endstep"
    
  end subroutine write_conf_frame

  subroutine read_conf_frame(st,fu, iostat)
    use geom
    use iso_fortran_env
    implicit none
    type(PcmcState) :: st
    integer, intent(in) :: fu
    integer,intent(out)::  iostat

    character(32) :: chardump
    character(1024) :: line
    integer :: i,k,Ndump,numero, fam, ntb
    real(8) :: dump,boxstar(3,3),vec3(3), rad
    iostat=0
    numero=0
    read(fu,'(A)', iostat=iostat) line
    if(iostat /= 0) then
      return
    end if
    read(line,'(A5,I8)') chardump, st%step
    read(fu,'(A9,I8)') chardump, ntb
    if(ntb>st%ntotbox) then
      call deallocate_boxes(st)
      st%ntotbox = ntb
      call allocate_boxes(st, st%Nbuff)
    end if
    read(fu,'(A)') chardump
    ! Nothing if not this token, for backwards compatibility
    if(chardump(1:13) == "#$config") then
      call read_changed_paras(st,fu,"frame")
      ! ends at a #$endconfig token
    end if
    st%ntotbox = ntb
    do k=1,st%ntotbox
      read(fu,'(A11)') chardump
      do i=1,3
          read(fu,'(3ES16.7)') st%boxes(k)%tens(i,1:3)
      end do
      dump = inverse_box(st%boxes(k)%tens,boxstar)
      
      st%boxes(k)%volume = dump
      call calc_metric(st%boxes(k)%tens, st%boxes(k)%met)
      
      read(fu,'(A21,I8)') chardump,Ndump
      if(Ndump>st%Nbuff) then
          write(error_unit,'(A,I0,A,I0,A)') 'ERROR: trajectory file: Please readjust the max number of particles to at least  '&
          ,Ndump, " (currently: ", st%Nbuff ,")."
          stop 7
      end if
      st%boxes(k)%n = Ndump
      do i=1,Ndump
          read(fu,'(A)') line
          read(line,'(3E16.8,I6,E16.8)') vec3(1:3),fam, rad
          st%boxes(k)%parts(i)%pos=matvec3_prod(boxstar,vec3)
          numero=numero+1
          st%boxes(k)%parts(i)%numero=numero
          if(fam<=0 .or. fam>st%dist%nfam) then
          write(error_unit,'(A,I0,A,I0,A)') "ERROR: trajectory file: Invalid family (" , fam, ") for Particule ", numero
          stop 7
          else
          st%boxes(k)%Parts(i)%famille = fam
          st%boxes(k)%Parts(i)%rayon = rad
          endif
      end do
      read(fu,*)
    end do
    ! May not be called if radii in the frame in the future
    ! call affect_radius_from_distrib(st)
    read(fu,'(A7)') chardump
  end subroutine read_conf_frame

  ! binary variant
  subroutine write_bconf_frame(st,fu)
    use iso_fortran_env
    implicit none
    type(PcmcState) :: st
    integer,intent(in) :: fu

    integer :: i,j,k
    integer(INT8) :: i2
    real(REAL32) :: vec(3),r
    real(8) :: bt(3,3)

    write(fu) bconfMag, 8, 0
    write(fu) st%step
    write(fu) st%ntotbox
    do k=1,st%ntotbox
      write(fu) bconfBox
      j=st%boxes(k)%n
      write(fu) j
      bt = st%boxes(k)%tens
      do i=1,3
        write(fu) bt(i,1:3)
      end do
      write(fu) bconfPos
      do i=1,j
        vec = real(st%boxes(k)%Parts(i)%pos, REAL32)
        write(fu) vec
      end do
      write(fu) bconfFam
      write(fu) st%dist%nfam
      if(st%dist%nfam <= 127) then
        do i=1,j
          i2 = int(st%boxes(k)%Parts(i)%famille, 1)
          write(fu) i2
        end do
      else
        do i=1,j
          write(fu) st%boxes(k)%Parts(i)%famille
        end do
      endif
      ! Write radii. Could be optional.
      write(fu) bconfRad
      do i=1,j
        r = real(st%boxes(k)%Parts(i)%rayon, REAL32)
        write(fu) r
      end do
      write(fu) bconfBoxEnd
    end do
  end subroutine

  subroutine read_bconf_frame(st,fu, iostat)
    use iso_fortran_env
    implicit none
    type(PcmcState) :: st
    integer,intent(in) :: fu
    integer,intent(out) :: iostat
   
    integer :: i,j,k,nfam
    integer :: hid,ri
    real(8) :: boxstar(3,3)
    real(REAL32) :: vec(3),r
    integer(INT8) :: i2
    read(fu, iostat=iostat) hid
    if(iostat /= 0) then
      return
    end if
    read(fu) ri, i
    if(hid /= bconfMag) call error("not a bconf file")
    read(fu) st%step
    read(fu) st%ntotbox
    do k=1,st%ntotbox
      read(fu) ri
      if(ri /= bconfBox) call error("unexpected token (expected box)")
      read(fu) j
      st%boxes(k)%n = j
      do i=1,3
        read(fu) st%boxes(k)%tens(i,1:3)
      end do
      st%boxes(k)%volume = inverse_box(st%boxes(k)%tens,boxstar)
      call calc_metric(st%boxes(k)%tens, st%boxes(k)%met)
      read(fu) ri
      if(ri /= bconfPos) call error("unexpected token (expected pos)")
      do i=1,j
        read(fu) vec
        st%boxes(k)%Parts(i)%pos = vec
      end do
      do while(.true.)
        read(fu) ri
        select case(ri)
        case(bconfFam)
          read(fu) nfam
          if(nfam <= 127) then
            do i=1,j
              read(fu) i2
              st%boxes(k)%Parts(i)%famille = i2
            end do
          else
            do i=1,j
              read(fu) st%boxes(k)%Parts(i)%famille
            end do
          end if
        case(bconfRad)
          do i=1,j
            read(fu) r
            st%boxes(k)%Parts(i)%rayon = r
          end do
        case(bconfBoxEnd)
          exit
        case default
          call error("unexpected token")
        end select
      end do
    end do
    ! May not be called if radii in the frame in the future
    ! call affect_radius_from_distrib(st)

  contains

    subroutine error(msg)
      character(*) :: msg
      write(error_unit,'(A,I0,A,A)') "Error, bconf file (unit=",fu,"): ",msg
      close(fu)
      stop 1
    end subroutine
  end subroutine

  subroutine write_changed_paras(st,ufp)
    type(PcmcState) :: st
    integer, intent(in) :: ufp
    write(ufp,'(A,G15.8)') "translation_interval = ",st%rpar%dx_trans*st%rpar%Lmin
    if(st%mode_box == 2) then ! to check
      write(ufp,'(A,G15.8)') "translation_interval_wall = ",st%rpar%dx_trans_wall*st%rpar%Lmin
    endif
  end subroutine

  subroutine read_changed_paras(st,ufp,fichpar)
    use readparam
    type(PcmcState) :: st
    integer, intent(in) :: ufp
    character(*), intent(in) :: fichpar
    character(len=:), allocatable :: key,val
    integer :: i,typ
    logical :: found
    type(ParamList) :: plist
    if(read_configuration(ufp, plist) < 0) then
      stop 1
    end if
    call start_iter_keyval(plist)
    do while(iter_keyval_in(plist))
      call get_current_keyval(plist, key, val)
      select case(key)
      case("translation_interval")
        read(val, *) st%rpar%dx_trans
      case("translation_interval_wall")
        read(val, *) st%rpar%dx_trans_wall
      end select
      call param_list_next(plist)
    end do
    !call configuration_finish(st)
  end subroutine read_changed_paras

end module state
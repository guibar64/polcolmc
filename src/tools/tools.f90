module tools
!! Various tools that can be called from the command line.
implicit none

contains

subroutine compute_gder(fichin, prefout,fich_dist, fich_par,dr_rdf, bstep, estep)
  !! Computes a radial distribution function from a trajectory.
  use state
  use correls
  use misc, only: short_initialisation
  use iso_fortran_env
  implicit none
  character(*), intent(in) :: fichin    !! trajectory file
  character(*), intent(in) :: prefout   !! prefix of output files
  character(*), intent(in) :: fich_dist !! distribution file
  character(*), intent(in) :: fich_par  !! configuration file
  real(8), intent(in) :: dr_rdf         !! space bin
  integer, intent(in) :: bstep, estep   !! range of steps to consider
  real(8) :: Lmin
  integer :: i,k, status, nframe, npframe, fu
  type(CorrelsRec),save :: res_rdf(max_boxes) ! rdf results and temps
  type(PcmcState) :: st
  type(Pmconf) :: fg

  call short_initialisation(st,fich_dist,fich_par)

  st%mode_box = get_mode_box(st%inp%smod_box)

  fg = open_pmconf(fichin, old=.true.)
  nframe=0
  npframe=0
  do
      call read_pmconf_frame(st,fg, .false., status)
      if(status == iostat_end) then
        exit
      else if(status /= 0) then
        write(error_unit, *) "Problem in input file. Exiting."
        stop 1
      end if
      nframe = nframe + 1

      if(nframe==1) then
        call allocate_multires(st)
        do k=1,st%ntotbox
            Lmin = st%boxes(k)%tens(1,1)
            do i=2, 3
              if(Lmin>st%boxes(k)%tens(i,i)) then
                  Lmin = st%boxes(k)%tens(i,i)
              end if
            end do
            call correls_init(res_rdf(k),st%dist%nfam,dr_rdf,0.D0,0.5D0*Lmin)

            call update_box_fam(st%boxes(k), st%dist)
            st%res(k)%popm(1:st%dist%nfam) = 0.d0
            st%res(k)%volm = 0.d0
          end do
      end if
      if(st%step > estep) exit
      if(st%step >= bstep) then
        npframe = npframe + 1
        do k=1,st%ntotbox
          call correls_new_step(res_rdf(k),st%boxes(k)%met,st%boxes(k)%n,st%boxes(k)%Parts,st%mode_box)
          st%res(k)%popm(1:st%dist%nfam) = st%res(k)%popm(1:st%dist%nfam) + &
            1.d0*st%boxes(k)%fam(1:st%dist%nfam)
          st%res(k)%volm = st%res(k)%volm + st%boxes(k)%volume
        end do
      end if
  end do
  call close_pmconf(fg)
  do k=1,st%ntotbox
      st%res(k)%volm = st%res(k)%volm/npframe
      st%res(k)%popm(1:st%dist%nfam) = st%res(k)%popm(1:st%dist%nfam)/npframe
      st%res(k)%Nm = sum(st%res(k)%popm(1:st%dist%nfam))
      call correls_end_step(res_rdf(k),st%res(k)%volm,st%res(k)%Popm)
  end do
  write(error_unit, '(A,I0,A)') "Processed ",npframe, " frames."
  do k=1,st%ntotbox
    call correls_reduit_g(res_rdf(k),st%res(k)%popm)
  end do

  open(newunit=fu,file=prefout)
  do k=1,st%ntotbox
      do i=1,res_rdf(k)%ng
        write(fu,'(ES16.7,1x,ES16.7)') (i-0.5D0)*dr_rdf, res_rdf(k)%gmoy(i)
      end do
      write(fu,'(A)') "&"
  enddo

  close(fu)


end subroutine

subroutine rescale_conf(fichin, fichout, factor,fich_dist,fich_par)
  !! Rescale a simulation box by `factor`.
  use state
  use misc, only: short_initialisation
  use iso_fortran_env
  implicit none
  character(*), intent(in) :: fichin    !! input positions
  character(*), intent(in) :: fichout   !! ouput positions
  character(*), intent(in) :: fich_dist !! distribution file
  character(*), intent(in) :: fich_par  !! configuration file
  integer :: k, status, nframe
  real(8) :: factor
  type(PcmcState) :: st
  type(Pmconf) :: fin, fout

  call short_initialisation(st,fich_dist,fich_par)

  fin = open_pmconf(fichin, old=.true.)
  nframe=0
  do
      call read_pmconf_frame(st,fin, .true., status)
      if(status == iostat_end) then
        exit
      else if(status /= 0) then
        write(error_unit, *) "Problem in input file. Exiting."
        stop 1
      end if
      nframe = nframe + 1
  end do
  call close_pmconf(fin)
  write(error_unit, '(A,I0,A)') "Read ",nframe, " frames."
  do k=1,st%ntotbox
    call box_set_tensor(st%boxes(k), st%boxes(k)%tens*(factor))
  end do
  fout = open_pmconf(fichout)
  call write_pmconf_frame(st,fout)
  call close_pmconf(fout)

end subroutine rescale_conf


subroutine confdouble(fichin, fichout, fich_dist, fich_par)
  !! Double a simulation box
  use state
  use misc, only: short_initialisation
  use iso_fortran_env
  implicit none
  character(*), intent(in) :: fichin    !! initial positions
  character(*), intent(in) :: fichout   !! final positions
  character(*), intent(in) :: fich_dist !! distribution file
  character(*), intent(in) :: fich_par  !! configuration file
  integer :: i,j,k, status, nframe
  real(8) :: alph(3)
  integer :: factor=2, px, py, pz, posk, ufout
  type(Box), allocatable :: mv2(:)
  type(PcmcState) :: st
  type(Pmconf) :: fin, fout

  call short_initialisation(st,fich_dist,fich_par)

  fin = open_pmconf(fichin, old=.true.)
  nframe=0
  do
      call read_pmconf_frame(st,fin, .true., status)
      if(status == iostat_end) then
        exit
      else if(status /= 0) then
        write(error_unit, *) "Problem in input file. Exiting."
        stop 1
      end if
      nframe = nframe + 1
  end do

  call close_pmconf(fin)
  write(error_unit, '(A,I0,A)') "Read ",nframe, " frames."

  allocate(mv2(st%ntotbox))
  do k=1,st%ntotbox
      allocate(mv2(k)%Parts(st%boxes(k)%n))
      mv2 = st%boxes
  end do

  call deallocate_boxes(st)
  do k=1, st%ntotbox
      st%boxes(k)%n = st%boxes(k)%n*factor**3
      allocate(st%boxes(k)%Parts(st%boxes(k)%n))
  end do


  do k=1,st%ntotbox
    call box_set_tensor(st%boxes(k), st%boxes(k)%tens*(factor))
    j = 1
    do px=0,factor-1 ; do py=0,factor-1 ; do pz=0,factor-1
      alph(1) = real(px,8)/real(factor)
      alph(2) = real(py,8)/real(factor)
      alph(3) = real(pz,8)/real(factor)
      do i=1,mv2(k)%n
          st%boxes(k)%Parts(j) = mv2(k)%Parts(i)
          do posk=1,3
            st%boxes(k)%Parts(j)%pos(posk) = mv2(k)%Parts(i)%pos(posk)/real(factor)+alph(posk)
          end do
          j = j + 1
      end do
    end do; end do ; end do
  end do
  st%dist%pop = st%dist%pop*factor**3
  st%step = 0 ! Because this is a new inst
  fout = open_pmconf(fichout)
  call write_pmconf_frame(st,fout)
  call close_pmconf(fout)

  open(newunit=ufout, file = trim(fich_dist)//".out")
  write(ufout,*) st%dist%nfam
  do i=1,st%dist%nfam
      write(ufout,'(I10,2ES16.7)') st%dist%pop(i), st%dist%rad(i), st%dist%rad(i)
  end do
  close(ufout)

end subroutine confdouble

subroutine do_mktables(dr,fichpar,fichdist)
  !! Create potential/force tables.
  use state
  use nrg
  implicit none
  real(8), intent(in) :: dr !! separation step
  character(*), intent(in) :: fichpar !! configuration file
  character(*), intent(in) :: fichdist !! distribution file
  integer :: stat
  type(PcmcState) :: st
  type(NrgRoutines) :: nr

  call read_cfg_file(st,fichpar)
  call init_distribution(st%dist, fichdist, stat)
  if(stat /= 0) return
  call init_modes(st)
  call nrg_init(st, nr)
  call make_xvg_tables(st, nr, dr)

contains

  subroutine make_xvg_tables(st,nr,dr_mktables)
    !! Make a bunch of tabulated potential files suitable for Gromacs
    type(PcmcState) :: st
    type(NrgRoutines) :: nr
    real(8), intent(in) :: dr_mktables
    character(256) :: fichmktables
    integer :: no,i,j,k,fu
    real(8) :: dump,dump1,rr

    no=floor((st%rc+1.D0)/dr_mktables)+3
    do i=1,st%dist%nfam
      do j=1,st%dist%nfam
          write(fichmktables,'("table_SPH",I0,"_SPH",I0,".xvg")') i,j
          open(newunit=fu,file=fichmktables)
          rr=0.D0
          write(fu,'(7ES16.7)') 0.D0,0.D0,0.D0,0.D0,0.D0,1.D7,1.D7
          do k=1,no
            rr=rr+dr_mktables
            if(rr>=st%rc) then
                dump=0.D0
                dump1=0.D0
            else
                dump=min(1.D7,nr%potential(rr,i,j,st))
                dump1=min(1.D7,nr%force(rr,i,j,st))
            end if
            write(fu,'(7ES16.7)') rr,0.D0,0.D0,0.D0,0.D0,dump,dump1
          end do
          close(fu)
      end do
    end do
  end subroutine

end subroutine

subroutine extract_conf(fichin, fichout, radrc, ipart,fich_dist, fich_par)
  !! Extract particles in a sphere around a particle and prints their positions in a XYZ file.
  use state
  use misc, only: dist2_mi, short_initialisation
  use iso_fortran_env
  implicit none
  character(*), intent(in) :: fichin    !! input positions
  character(*), intent(in) :: fichout   !! extracted positions
  character(*), intent(in) :: fich_dist !! distribution file
  character(*), intent(in) :: fich_par  !! configuration file
  real(8), intent(in) :: radrc          !! radius of the sphere
  integer, intent(in) :: ipart          !! index of the particle
  character(*), parameter :: exename="polcolmc"
  integer :: i,j,k, status, nframe, is=0
  real(8) :: relpos(3), rr, radrc2
  type(Particule), allocatable :: selection(:)
  type(PcmcState) :: st
  type(Pmconf) :: fin
  integer :: ufout

  call short_initialisation(st,fich_dist,fich_par)

  fin = open_pmconf(fichin, old=.true.)
  nframe=0
  do
      call read_pmconf_frame(st,fin, .true., status)
      if(status == iostat_end) then
        exit
      else if(status /= 0) then
        write(error_unit, *) "Problem in input file. Exiting."
        stop 1
      end if
      nframe = nframe + 1
  end do
  call close_pmconf(fin)
  write(error_unit, '(A,I0,A)') "Read ",nframe, " frames."

  st%mode_box = get_mode_box(st%inp%smod_box)

  radrc2 = radrc*radrc
  do k=1,st%ntotbox
      do i=1,st%boxes(k)%n
        if(st%boxes(k)%parts(i)%numero == ipart) then
            allocate(selection(st%boxes(k)%n))
            is = 1
            selection(is) = st%boxes(k)%parts(i)
            do j=1,st%boxes(k)%n
              if(j==i) cycle
              rr = dist2_mi(selection(1)%pos-st%boxes(k)%parts(j)%pos, relpos, st%boxes(k)%met,st%mode_box)
              if(rr <= radrc2) then
                  is = is + 1
                  selection(is) = st%boxes(k)%parts(j)
                  selection(is)%pos = matvec3_prod(st%boxes(k)%tens,relpos)
              end if
            end do
            exit
        end if
      end do
  end do
  write(error_unit, '(I0," particles selected around particle ",I0," (",I0,") in a st%dist%rad of ", G15.8)') &
        is, ipart, selection(1)%famille, radrc
  if(is > 0) selection(1)%pos = 0.d0
  open(newunit=ufout,file=fichout)
  call write_xyz(ufout, is, selection, st%step)
  close(ufout)

contains

  subroutine write_xyz(uf, np, partv, step)
    implicit none
    integer :: uf, step, np
    type(Particule) :: partv(np)
    intent(in) :: uf, np, partv
    integer :: i
    real(8) :: pos(3)
    character(6) :: chardump
    write(uf, *) np
    write(uf,'(A,A,A,I0)') "Generated by ", exename,", step = ",step
    do i=1,np
        pos = partv(i)%pos
        write(chardump, '("SPH",I0)') partv(i)%Famille
        write(uf,'(A6,3E16.8)') chardump , pos
    end do
  end subroutine write_xyz

end subroutine


subroutine convert_conf(nin, fichin, fichout, fich_dist, fich_par, iformat, oformat, ibox)
  !! Convert a trajectory file to another format. Works on **one single** box.
  use state
  use trr
  use misc, only: short_initialisation
  use iso_fortran_env
  implicit none
  integer, intent(in) :: nin             !! number of input files 
  character(*), intent(in) :: fichin(nin) !! input trajector-y-ies
  character(*), intent(in) :: fichout    !! output trajectory
  character(*), intent(in) :: fich_dist  !! distribution file
  character(*), intent(in) :: fich_par   !! configuration file
  character(*), intent(in) :: iformat    !! input format
  character(*), intent(in) :: oformat    !! output format
  integer, intent(in) :: ibox            !! box index 

  integer, parameter :: FCONF=1,FBCONF=2,FXYZ=3,FGRO=4,FTRR=5
  character(*), parameter :: exename="polcolmc"
  integer :: status, nframe, ifmt, ofmt, ntotbox,k
  type(PcmcState) :: st
  type(Pmconf) :: fin, fout
  integer, allocatable :: ufin(:)
  integer :: ufout
  
  if(nin<1) then
     write(error_unit,'(A)') "Error: no input trajectory. Nothing to do."
    stop 0
  end if

  if(len_trim(iformat) == 0) then
    ifmt = guess_format(fichin(1)) ! First file determine format
  else
    ifmt = select_format(iformat)
  end if
  if(len_trim(oformat) == 0) then
    ofmt = guess_format(fichout)  
  else
    ofmt = select_format(oformat)
  end if
  
  if(ifmt == FTRR) then
    ! trr files does not store information about type, could we assume the particles
    ! are in the 'right' order ?
    write(error_unit,'(A)') "Error: trajectory formats 'trr' not available for reading."
    stop 0
  end if

  call short_initialisation(st,fich_dist,fich_par)

  select case(ifmt)
  case(FCONF, FBCONF)
    ntotbox = st%ntotbox ! ⚠ The file could claim otherwise
  case(FXYZ, FGRO, FTRR)
    ntotbox = min(nin, max_boxes)
  end select
  if(ntotbox<nin) write(error_unit,'(A,I0,A)') "Note: using only ",ntotbox," files"

  if(ntotbox<ibox) then
    write(error_unit,'(A,I0)') "Error: Box index out of range: ",ibox
    stop 1
  end if
  allocate(ufin(ntotbox))

  select case(ifmt)
  case(FCONF, FBCONF)
    ! These formats support multiple boxes so other input files are ignored.
    fin = open_pmconf(fichin(1), old=.true.)
  case(FXYZ, FGRO)
    ! no box info for these formats
    if(st%ntotbox < ntotbox) then
      call deallocate_boxes(st)
      call allocate_boxes(st, st%Nbuff) ! Could try to preread the files instead
    end if
    call init_box_from_configuration(st,st%inp%boxLx,st%inp%boxLy,st%inp%boxLz,&
      st%inp%boxalpha,st%inp%boxbeta,st%inp%boxgamma,st%inp%phi_inp)
    do k=1,ntotbox
      open(newunit=ufin(k), file = fichin(k), status='old', action='read')
    end do
  case(FTRR)
    do k=1,ntotbox
      call open_trr(ufin(k), fichin(k))
    end do
  end select

  select case(ofmt)
  case(FCONF)
    fout = open_pmconf(fichout)
  case(FBCONF)
    fout = open_pmconf(fichout, bin=.true.)
  case(FXYZ, FGRO)
    open(newunit=ufout, file=fichout, action='write')
  case(FTRR)
    call open_trr(ufout, fichout, w=.true.)
  end select

  nframe=0
  do
    select case(ifmt)
    case(FCONF, FBCONF)
      call read_pmconf_frame(st, fin, .false., status) 
      ntotbox = st%ntotbox
    case(FXYZ)
      do k=1,ntotbox
        call read_xyz(ufin(k), st%boxes(k), st%step, status)
      end do
    case(FTRR)
      do k=1,ntotbox
        call read_trr(ufin(k), st%boxes(k), st%step, status)
      end do
    case(FGRO)
      do k=1,ntotbox
        call read_gro(ufin(k), st%boxes(k), st%step,status)
      end do
    end select
    if(status == iostat_end) then
      exit
    else if(status /= 0) then
      write(error_unit, '(A,I0,A)') "Error: Problem in input file (code ",status,")."
      stop 1
    end if

    ! Type may change from frame to frame    
    select case(ifmt)
    case(FCONF,FBCONF)
      ! Radii are included inf these formats
    case(FXYZ,FTRR,FGRO)
      call affect_radius_from_distrib(st)
    end select

    nframe = nframe + 1

    select case(ofmt)
    case(FCONF,FBCONF)
      call write_pmconf_frame(st,fout,no_cfg=.true.)
    case(FXYZ)
      call write_xyz(ufout, st%boxes(ibox),st%step)
    case(FGRO)
      call write_gro(ufout, st%boxes(ibox),st%step)
    case(FTRR)
      call write_trr(ufout, st%boxes(ibox),st%step)
    end select
  end do
  
  select case(ifmt)
    case(FCONF, FBCONF)
    call close_pmconf(fin)
  case default
    do k=1,ntotbox
      close(ufin(k))
    end do
  end select
  select case(ofmt)
    case(FCONF, FBCONF)
    call close_pmconf(fout)
  case default
    close(ufout)
  end select

  write(error_unit, '(A,I0,A)') "Processed ",nframe, " frames."

  contains

  integer function guess_format(fileName)
    character(*) :: fileName
    integer :: iend, start
    do iend=len(fileName),1,-1
      if(fileName(iend:iend) /= ' ') exit
    end do
    if(iend == 1) return
    do start=iend,1,-1
      if(fileName(start:start) == '.') exit
    end do
    start=start+1
    if(iend-start<0) return
    guess_format = select_format(fileName(start:iend))
  end function

  integer function select_format(format)
    character(*), intent(in) :: format
    select case(trim(format))
      case("xyz")
        select_format = FXYZ
      case("gro")
        select_format = FGRO
      case("trr")
        select_format = FTRR
      case("inst-bin")
        select_format = FBCONF
      case default
        select_format = FCONF
    end select
  end function

  subroutine write_xyz(uf, b,step)
    implicit none
    integer :: uf,step
    type(Box) :: b
    intent(in) :: uf, b
    integer :: i
    real(8) :: pos(3)
    character(6) :: chardump
    write(uf, *) b%n
    write(uf,'("L = ",ES16.7,", step = ",I10)') b%tens(1,1),step
    do i=1,b%n
        pos = matvec3_prod(b%tens,  b%parts(i)%pos)
        write(chardump, '("SIL",I0)') b%parts(i)%Famille
        write(uf,'(A6,3ES16.7)') chardump , pos
    end do
  end subroutine write_xyz

  subroutine read_xyz(uf, b, step,stat)
    implicit none
    integer :: uf, stat,step
    type(Box) :: b
    intent(in) :: uf
    integer :: i, statdesc
    real(8) :: pos(3), ibv(3,3)
    character(6) :: chardump
    character(32) :: stemp1,stemp2
    character(3) prefname
    character(256) :: info
    read(uf, *, iostat=stat) b%n
    if(stat /= 0) return
    ! cubic box for now
    read(uf, '(A)') info
    if(len(info)> 0 .and. info(1:1) == 'L') then
      ! Can fail
      read(info,'(A4,ES16.7,A9,I10)', iostat=statdesc) stemp1,b%tens(1,1),stemp2,step
      if(statdesc /= 0) write(error_unit, '(A,A,A)') "Warning: can not extract box length from '", trim(info),"'"
    end if
    b%tens(2,2) =  b%tens(1,1)
    b%tens(3,3) =  b%tens(1,1)
    pos(1) = inverse_box(b%tens, ibv)
    do i=1,b%n
       read(uf,'(A6,3ES16.7)') chardump , pos
       read(chardump, '(A3,I3)') prefname, b%parts(i)%Famille
       b%parts(i)%pos = matvec3_prod(ibv, pos)
    end do
  end subroutine read_xyz

  subroutine write_gro(uf, b,step)
    implicit none
    integer :: uf,step
    type(Box) :: b
    intent(in) :: uf, b
    integer :: i
    real(8) :: pos(3)
    character(7) :: chardump
    write(uf, '(I0)') b%n
    do i=1,b%n
        pos = matvec3_prod(b%tens,  b%parts(i)%pos)
        write(chardump,'("SIL",I0)') b%parts(i)%Famille
        write(uf, '(I5,A7,A2,I5,1x,F7.3,1x,F7.3,1x,F7.3,1x,F7.3,1x,F7.3,1x,F7.3)') i, chardump, "Si", i,pos, 0.0,0.0,0.0
    end do
    write(uf, '(1x,F9.5,1x,F9.5,1x,F9.5)') b%tens(1,1),  b%tens(2,2),  b%tens(3,3)
  end subroutine write_gro

  subroutine read_gro(uf, b, step,stat)
    implicit none
    integer :: uf, stat,step
    type(Box) :: b
    intent(in) :: uf
    integer :: i,j
    real(8) :: pos(3), ibv(3,3), buffer
    character(6) :: chardump
    character(3) prefname
    character(256) :: info
    character(1024) line
    character(2) :: res
    allocatable :: buffer(:,:)
    read(uf, *, iostat=stat) b%n
    if(stat /= 0) return
    allocate(buffer(3,b%n))
    read(uf,'(A)') info ; call parse_info(info,step)
    do i=1,b%n
      read(uf,'(A)') line
      read(line,'(I5,A7,A2,I5,F8.3,F8.3,F8.3)') j,chardump , res ,j,buffer(i,1), buffer(i,2), buffer(i,3)
      read(chardump, '(A3,I3)') prefname, b%parts(i)%Famille
      b%parts(i)%pos = matvec3_prod(ibv, pos)
    end do
    b%tens = 0  ! And if not cubic ?
    read(uf, '(F10.5,F10.5,F10.5)') b%tens(1,1),  b%tens(2,2),  b%tens(3,3)
    pos(1) = inverse_box(b%tens, ibv)
    do i=1,b%n
      b%parts(i)%pos = matvec3_prod(ibv, buffer(1:3,i))
    end do
  end subroutine read_gro


  subroutine parse_info(line, step)
    character(*), intent(in) :: line
    integer, intent(out) :: step
    integer :: spos
    character(128) :: val
    step = 0
    spos = find_keyval(line,"step",val,1)
    if(spos > 0) then
      read(val,*) step
    end if
  end subroutine

  integer function find_keyval(line, key, val,start) result(pos)
    character(*), intent(in) :: line,key
    character(*), intent(out) :: val
    integer, intent(in) :: start
    integer :: i,j,k,l,ival
    if(len(key)<1) then
      pos = 0
      return
    end if
    do i=start,len(line)-len(key)
      if(line(i:i+len(key)-1) == key ) exit
    end do
    if(i >= len(line)-len(key)) then
      pos = 0
      return
    end if
    i = i + len(key)
    do j=i,len(line)
      if(line(j:j) == '=') then
        exit
      else if(line(j:j) /= ' ') then ! could use other whitespaces
        pos = 0
        return
      end if
    end do
    if(j>=len(line)) then
      pos = 0
      return
    end if
    do k = j+1,len(line)
      if(line(k:k) /= ' ') exit
    end do
    ival = 1
    do l=k,len(line)
      if(line(l:l) == ' ') exit
      val(ival:ival) = line(l:l)
      ival = ival + 1
    end do
    val(ival:) = ' '
    pos = l
  end function

end subroutine

subroutine analyze_bondorder(fichin, fich_dist, fich_par, ibox, fich_bocfg)
  !! Performs bond-order analysis on a trajectory.
  use state
  use bondorder
  use misc, only: short_initialisation, dist2_mi
  use iso_fortran_env
  use readparam
  implicit none
  character(*), intent(in) :: fichin     !! input trajectory
  character(*), intent(in) :: fich_dist  !! distribution file
  character(*), intent(in) :: fich_par   !! trajectory file
  character(*), intent(in), optional :: fich_bocfg !! bond-order configuration
  integer, intent(in) :: ibox            !! box index
  real(8) :: Lmin
  integer :: i,k, s, status, nframe, ufc
  type(PcmcState) :: st
  type(Pmconf) :: fg
  type(ParamList) :: config
  type(OrderData), allocatable :: x(:)
  integer :: np
  type(BondResults) :: results
  real(8) :: rpos(3), pos(3), r2
  real(8), allocatable :: inip(:,:), relp(:,:)
  
  call short_initialisation(st,fich_dist,fich_par)

  st%mode_box = get_mode_box(st%inp%smod_box)
  np = sum(st%dist%pop(:))
  allocate(x(np), relp(3,np), inip(3,np))
  if(present(fich_bocfg)) then
    open(newunit=ufc, file=fich_bocfg, status='old')
    if(read_configuration(ufc, config, "bondorder") /= 0) then
      write(error_unit, '(A)') "Error in config file '", trim(fich_bocfg), "'."
      stop 1
    end if

    call bondorder_initialize(np, x, results, st%dist, config)
  else
    call bondorder_initialize(np, x, results, st%dist, st%config)
  end if
  fg = open_pmconf(fichin, old=.true.)
  nframe=0
  do s=1,10000 ! TODO: multiple frame analysis
    call read_pmconf_frame(st,fg, .false., status)
    if(status == iostat_end) then
      exit
    else if(status /= 0) then
      write(error_unit, *) "Problem in input file. Exiting."
      stop 1
    end if
    nframe = nframe + 1

    if(nframe==1) then
      call allocate_multires(st)
      do k=1,st%ntotbox
        Lmin = st%boxes(k)%tens(1,1)
        do i=2, 3
          if(Lmin>st%boxes(k)%tens(i,i)) then
            Lmin = st%boxes(k)%tens(i,i)
          end if
        end do

        call update_box_fam(st%boxes(k), st%dist)
      end do
    else
    end if
    call bondorder_do_step(st%boxes(ibox)%n, x, results, st%boxes(ibox), st%step, st%dist)
  end do
  call close_pmconf(fg)
  call bondorder_finalize(np, x, results)
end subroutine

subroutine compute_msd(fichin, fichout,fich_dist, fich_par, ibox, bstep, estep)
  !! Computes the mean square displacement from a trajectory.
  use state
  use correls
  use misc, only: short_initialisation, dist2_mi
  use readparam, only: get_param
  use iso_fortran_env
  implicit none
  character(*), intent(in) :: fichin    !! trajectory file
  character(*), intent(in) :: fichout   !! output file
  character(*), intent(in) :: fich_dist !! distribution file
  character(*), intent(in) :: fich_par  !! configuration file
  integer, intent(in) :: ibox           !! box index
  integer, intent(in) :: bstep          !! first step to consider
  integer, intent(in) :: estep          !! last step to consider

  integer :: nsamples, psample, isample, lmax, nl
  integer :: i,k,l, status, nframe, step_init, fam
  integer :: ufo=47, uft = 49
  type(CorrelsRec),save :: res_rdf(max_boxes) ! rdf results and temps
  type(PcmcState) :: st
  type(Pmconf) :: fg
  real(8) :: dpos(3), dposout(3),r2, msdm
  real(8), allocatable :: x(:,:,:), msd(:), popm(:)
  integer, allocatable :: steps(:)
  type MetContainer
    real(8) :: data(3,3)
  end type
  type(MetContainer), allocatable :: mets(:)
  ! for diffusion
  real(8) :: phi,x2m,xm,ym,xym,norm,Dfake, vpm, volm
  integer :: stepc
  character(32) :: fmt_msd
  character(len=:), allocatable :: val
  logical :: opened

  call short_initialisation(st,fich_dist,fich_par)
  st%mode_box = get_mode_box(st%inp%smod_box)

  if(get_param(st%config, "msd.number_of_samples", val)) then
    read(val,*) nsamples
  else
    nsamples = 1
  end if
  if(get_param(st%config, "msd.sample_period", val)) then
    read(val,*) psample
  else
    psample = 1
  end if


  fg = open_pmconf(fichin, old=.true.)
  open(newunit=ufo, file=fichout)

  allocate(msd(st%dist%nfam),popm(st%dist%nfam))
  write(fmt_msd,'(A,I0,A)') "(I10,",st%dist%nfam+1,"ES16.7)"
  
  xm=0
  ym=0
  x2m=0
  xym=0
  norm=0
  volm = 0.
  vpm = 0.
  nframe=0
  l = 1
  step_init = 0
  ! Det number of frames
  do
    call read_pmconf_frame(st,fg, .false., status)

    if(status == iostat_end .or. st%step > estep) then
      exit
    else if(status /= 0) then
      write(error_unit, *) "Problem in input file. Exiting."
      stop 1
    end if
    if(nframe == 0) step_init = st%step
    if(st%step >= bstep .and. st%step <= estep) then
      nframe = nframe + 1
    end if
  end do
  nframe = nframe - 1
  rewind fg%unit
  call close_pmconf(fg)
  ! 🤔 What if number of parts changes ?
  allocate(x(3, ubound(st%boxes(ibox)%parts,1), 0:nframe), mets(0:nframe))
  allocate(steps(0:nframe))
  
  fg = open_pmconf(fichin, old=.true.)
  popm = 0.0_8
  status = 0 ! why do I have to reset that variable ?
  ! Read entire trajectory
  k = 0
  do
    call read_pmconf_frame(st,fg, .false., status)

    if(status == iostat_end .or. st%step > estep) then
      exit
    else if(status /= 0) then
      inquire(unit=fg%unit, opened=opened)
      write(error_unit, *) "Problem in input file, frame ",k, ". Exiting."
      stop 1
    end if
    if(st%step >= bstep .and. st%step <= estep) then
      do i=1,st%boxes(ibox)%n
        x(1:3, i, k) = st%boxes(ibox)%parts(i)%pos
        fam = st%boxes(ibox)%parts(i)%famille
        popm(fam) = popm(fam) + 1._8
        vpm = vpm + st%boxes(ibox)%parts(i)%rayon**3
      end do
      if(k == 0) step_init = st%step
      steps(k) = st%step - step_init
      mets(k)%data = st%boxes(ibox)%met
      volm = st%boxes(ibox)%volume + volm 

      k = k + 1
    end if
  end do
  popm = popm / (nframe+1)
  vpm = vpm / (nframe+1)
  volm = volm / (nframe+1)
  do k=1, nframe
    msd = 0.0_8
    msdm = 0.0_8
    lmax = min((nsamples-1)*psample,nframe-k)
    nl = 0
    do l=0,lmax,psample
      nl = nl + 1
      do i=1,ubound(x,2)
        dpos = x(1:3, i, k+l) - x(1:3, i,l)
        r2 = dist2_mi(dpos, dposout, mets(k)%data, st%mode_box)
        fam = st%boxes(ibox)%parts(i)%famille
        msd(fam) = msd(fam) + r2
        msdm = msdm + r2
      end do
    end do
    msd = msd / (popm*nl)
    msdm = msdm / (ubound(x,2)*nl)
    stepc = steps(k)
    ym=ym+msdm*stepc
    xym=xym+stepc*msdm*stepc
    xm=xm+stepc*stepc
    x2m=x2m+(1.D0*stepc)**3
    norm=norm+stepc

    write(ufo,fmt=fmt_msd) stepc,msdm,msd(:)
  end do
  call close_pmconf(fg)
  close(ufo)
  write(error_unit, '(A,I0,A)') "Processed ",nframe+1, " frames."

  ym=ym/norm
  xym=xym/norm
  xm=xm/norm
  x2m=x2m/norm
  Dfake=(xym-xm*ym)/(x2m-xm*xm)

  phi = 4._8*PI/3._8 * vpm / volm

  write(output_unit, '(A)') "phi           Dfake"
  write(output_unit,'(2ES16.7)') phi,Dfake

 
end subroutine

subroutine prune_conf(fich_dist,fich_par, fichin, fichout, period, last_step)
  !! Reduce a trajectory to every `period` steps.
  use state
  use misc, only: short_initialisation
  use iso_fortran_env
  implicit none
  character(*), intent(in) :: fichin    !! input positions
  character(*), intent(in) :: fichout   !! ouput positions
  character(*), intent(in) :: fich_dist !! distribution file
  character(*), intent(in) :: fich_par  !! configuration file
  integer, intent(in) :: period, last_step

  integer :: k, status, nframe, last, step_init, nwritten
  type(PcmcState) :: st
  type(Pmconf) :: fin, fout
  if(last_step < 0) then 
    last = huge(1)
  else
    last = last_step
  end if
  call short_initialisation(st,fich_dist,fich_par)

  fin = open_pmconf(fichin, old=.true.)
  fout = open_pmconf(fichout, bin=.true.)
  step_init = 0
  nframe=0
  nwritten = 0
  do
      call read_pmconf_frame(st,fin, .true., status)
      if(status == iostat_end) then
        exit
      else if(status /= 0) then
        write(error_unit, *) "Problem in input file. Exiting."
        stop 1
      end if
      if(nframe == 0) step_init = st%step
      if(st%step > last) exit
      nframe = nframe + 1
      if(modulo(st%step - step_init, period) == 0) then
        call write_pmconf_frame(st, fout)
        nwritten = nwritten + 1
      end if
  end do
  write(error_unit, '(A,I0,A)') "Read ", nframe, " frames."
  write(error_unit, '(A,I0,A)') "Written ", nwritten, " frames."
  call close_pmconf(fin)
  call close_pmconf(fout)

end subroutine


end module tools

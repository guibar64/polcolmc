module readparam
  !! Reading configuration files composed of  `key=value`lines.
  !! 
  !!## Example
  !!
  !!```ini
  !!key = value
  !!
  !![section]
  !!key = value
  !!```
  !!
  !!```fortran
  !!use readparam
  !!type(ParamList) :: list
  !!integer :: fconf
  !!character(len=:), allocatable :: value
  !!open(newunit=fconf, file="myconf.cfg")
  !!call read_configuration(fconf, list)
  !!if(get_param(list, "key", value)) then
  !!  write(*,*) "'key' found, value = ", value
  !!else
  !!  write(*,*) "key not found"
  !!end if
  !!if(get_param(list, "section.key", value)) then
  !!  write(*,*) "in section 'section': 'key' found, value = ", value
  !!else
  !!  write(*,*) "key not found"
  !!end if
  !!```
  use iso_fortran_env
  implicit none
  
  integer, parameter, private :: param_list_min=128, param_list_chunk=100

  type Param
    character(len=:), allocatable :: key, val
    logical :: accessed
  end type Param

  type ParamList
    integer :: len
    integer :: idx
    type(Param), allocatable :: pars(:)
  end type ParamList

  private :: Param, add_par

contains

subroutine init_param_list(list)
  type(ParamList), intent(out) :: list
  list%len = 0
  list%idx = 0
  allocate(list%pars(param_list_min))
end subroutine init_param_list

subroutine clear_param_list(list)
  type(ParamList), intent(inout) :: list
  list%len = 0
  list%idx = 0
  deallocate(list%pars)
end subroutine

subroutine start_iter_keyval(list)
  type(ParamList), intent(inout) :: list
  list%idx = 1
end subroutine start_iter_keyval

logical function iter_keyval_in(list)
  type(ParamList), intent(in) :: list
  iter_keyval_in = list%idx <= list%len
end function iter_keyval_in

subroutine param_list_next(list)
  type(ParamList) :: list
  list%idx = list%idx + 1
end subroutine param_list_next

subroutine get_current_keyval(list, key,val)
  type(ParamList), intent(in) :: list
  character(len=:), allocatable :: key,val
  integer :: kl,vl
  key = list%pars(list%idx)%key
  val = list%pars(list%idx)%val
end subroutine

subroutine access_current(list)
  type(ParamList), intent(inout) :: list
  list%pars(list%idx)%accessed = .true.
end subroutine access_current

subroutine add_par(list, key, val)
  type(ParamList), intent(inout) :: list
  character(*), intent(in) :: key,val
  integer :: nlen, ocap
  type(Param), allocatable :: temp(:)
  nlen = list%len
  nlen = nlen + 1
  ocap = ubound(list%pars,1)
  if(nlen > ocap) then
    call move_alloc(list%pars, temp)
    allocate(list%pars(ocap+param_list_chunk))
    list%pars(1:ocap) = temp(1:ocap)
    deallocate(temp)
  end if
  list%pars(nlen)%key = key
  list%pars(nlen)%val = val
  list%pars(nlen)%accessed = .false.
  list%len = nlen
end subroutine add_par

integer function read_configuration(unit, list, initial_section)
  !! reads a cfg file
  type(ParamList), intent(inout) :: list !! parameter list
  integer, intent(in) :: unit  !! file unit to read from
  character(*), intent(in), optional :: initial_section !! sets the initial section for otherwise sectionless keys
  character(512) :: key,val
  integer :: i,kl,vl,typ, statut, c,j, k, keyLen, valLen, nl, fp
  character(1024) :: line
  character(256) :: section
  integer, parameter :: maxKeyLen = 256, maxValLen = 512
  logical :: insec
  read_configuration = 0
  do i = 1, len(val)
    val(i:i) = ' '
  end do
  do i=1,len(line)
    line(i:i) = ' '
  end do
  call init_param_list(list)
  if(present(initial_section)) then
    section = initial_section
    insec = .true.
  else
    insec = .false.
  end if
  do while(.true.)
    read(unit,'(a)',iostat=statut) line
    nl = nl + 1
    if(statut== iostat_end) then
      exit
    else if(statut /= 0) then
      write(error_unit, '(A,I0)') "Error: Expected line ", nl
      read_configuration = -1
      exit
    end if
    do fp=1,len(line)
      if(line(fp:fp) /= ' ') exit
    end do
    if(fp>len(line)) then
      ! empty 
    else if(line(fp:fp)=="!" .or. line(fp:fp)=="#") then
      if(line(fp+1:fp+10)=='$endconfig') then
        exit
      else
        ! comment
      end if
    else if(line(fp:fp) == '[') then
      ! section
      do i=fp+1, len(line)
        if(line(i:i) /= ' ') exit
      end do
      do c=i,len(line)
        if(line(c:c) == ']') exit
      end do
      if(c>len(line)) then
        write(error_unit, '(A,I0,A)') "Error line ", nl,": Expected ']'"
      end if
      do j=c,len(line)
        if(line(j:j) /= ' ') exit
      end do
      section = line(i:j-1)
      insec = .true.
    else if(len_trim(line) == 0) then
      ! Empty line
    else
      do i=fp,fp+maxKeyLen-1
        if(line(i:i)=="=") exit
      end do
      if(i>=fp+maxKeyLen) then
        write(error_unit, '(A,I0,A,I0,A)') "Error line ", nl, ": Expected '=' line or key too longlength(>",maxKeyLen," chars)"
        read_configuration = -1
        exit
      end if
      do c=i-1,fp
        if(line(c:c) /= ' ') exit
      end do
      keyLen = min(c-fp+1,len(key))
      if(insec) then
        key = trim(section) // "." // line(fp:fp+keyLen-1)
        keyLen = keyLen + len_trim(section) + 1
      else
        key = line(fp:fp+keyLen-1)
      end if
      i = i + 1
      do j=i,maxValLen
        if(line(j:j) /= ' ') exit
      end do
      if(j>maxValLen) then
        write(error_unit, '(A,I0,A,I0,A)') "Error line ", nl,": value too long (>", maxValLen, " chars) after '='"
        read_configuration = -1
        exit
      end if
      if(line(j:j) == '"') then
        j = j + 1
        do k=j, len(line)
          if(line(k:k) == '"') exit
        end do
        if(k>len(line)) then
          write(error_unit, '(A,I0,A)') "Error line ", nl,': expected ''"'' before end of line.'
          read_configuration = -1
          exit
        end if
        k = k - 1 
      else
        do k=j, len(line)
          if(line(k:k) == ' ') exit
        end do
        k = k - 1
      end if
      valLen = min(k - j + 1, maxValLen)
      val = line(j:k)
      call add_par(list, key(1:keyLen), val(1:valLen))
    end if
  end do
end function

logical function get_param(list, key, val)
  type(ParamList) :: list
  character(*), intent(in) :: key
  character(len=:), allocatable :: val
  integer :: i
  get_param = .false.
  do i=list%len, 1, -1
    if(list%pars(i)%key(1:len_trim(list%pars(i)%key)) == key) then
      list%pars(i)%accessed = .true.
      val = list%pars(i)%val
      get_param = .true.
      return
    end if
  end do
end function

subroutine warn_about_unkown_params(unit,list)
  integer, intent(in) :: unit
  type(ParamList), intent(in) :: list
  integer :: i
  do i=1, list%len
    if(.not.list%pars(i)%accessed) then
      write(unit,'(A,A,A)') "Warning: Unkown parameter '", list%pars(i)%key, "'."
    end if
  end do
end subroutine

function read_logical(val) result(res)
  logical :: res
  character(*), intent(in) :: val
  integer :: lv
  lv = len_trim(val)
  select case(val(1:lv))
  case('y',"yes","true","on","oui")
    res = .true.
  case default
    res = .false.
  end select
end function read_logical

end module readparam

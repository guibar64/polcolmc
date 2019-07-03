module trr
!! Reading and writing gromacs TRR trajectories.
use iso_fortran_env
implicit none

contains

#ifdef NO_CONVERSPEC_SUPPORT
subroutine open_trr(fu, file, w)
  integer, intent(out) :: fu
  character(*) :: file
  logical, intent(in), optional :: w
  fu = 0
  write(error_unit,'(A)') "Error: This program was not built with support for trr."
  stop 1
end subroutine open_trr
#else
subroutine open_trr(fu, file, w)
  integer, intent(out) :: fu
  character(*), intent(in) :: file
  logical, intent(in), optional :: w
  ! non standard convert specifier supported by gfortran and ifort (at least)
  if(present(w) .and. w) then
    open(newunit=fu,file=file,status='replace',access="stream",form="unformatted", convert="big_endian", action='write')
  else
    open(newunit=fu,file=file,status='old', access="stream",form="unformatted", convert="big_endian", action='read')
  end if
end subroutine
#endif

subroutine write_trr(uf, b,step)
  !! write box info and positions to an opened file in trr format.
  use state, only: box
  use geom, only: matvec3_prod
  implicit none
  integer, intent(in) :: uf !! file descriptor
  integer, intent(in) :: step !! current step
  type(Box), intent(in) :: b  !! simulation box
  integer :: i,j
  real(8) :: pos(3)
  integer(4) :: junk(7)
  junk = 0
  write(uf) 1993
  write(uf) 13
  write(uf) 12
  write(uf) "GMX_trn_file"
  write(uf) 0
  write(uf) 0
  write(uf) 9*8  ! box size
  write(uf) 0
  write(uf) 0
  write(uf) 0
  write(uf) 0
  write(uf) b%n*8*3 ! xsize
  write(uf) 0  !vsize
  write(uf) 0  ! fsize
  write(uf) b%n
  write(uf) step
  write(uf) 0
  write(uf) real(step, 8)
  write(uf) 0.d0  ! lambda
  do i=1,3 ; do j=1,3
    write(uf) b%tens(j,i)
  end do ;enddo
  do i=1,b%n
      pos = matvec3_prod(b%tens,  b%parts(i)%pos)
      do j=1,3
        write(uf) pos(j)
      end do
  end do
end subroutine write_trr

subroutine read_trr(uf, b,step, status)
  !! write box info and positions to an opened file in trr format.
  use state, only: box
  use geom, only: matvec3_prod, inverse_box, calc_metric
  implicit none
  integer :: uf   !! file descriptor
  integer :: step !! current step
  integer :: status !! set to /= 0 on error
  type(Box) :: b    !! simulation box
  intent(in) :: uf
  intent(out) :: status, step
  integer :: i,j
  real(8) :: pos(3), ld, bt(3,3)
  real(8) :: boxstar(3,3)
  character(12) :: chardump
  integer(4) :: junk(7), mag, v1, v2, box_size, xsize,vsize,fsize,np
  status = 0
  junk = 0
  read(uf, iostat=status) mag
  if(status /= 0) then
    return
  end if
  if(mag /= 1993) then
    write(*,'(A)') "Invalid trr file"
    status = 1
    return
  end if
  read(uf) v1
  read(uf) v2
  read(uf) chardump
  read(uf) i
  read(uf) i
  read(uf) box_size  ! box size
  read(uf) i
  read(uf) i
  read(uf) i
  read(uf) i
  read(uf) xsize
  read(uf) vsize
  read(uf) fsize
  read(uf) np
  read(uf) step
  read(uf) i
  read(uf) ld  ! time
  read(uf) ld  ! lambda
  if(box_size == 9*8) then
    do i=1,3 ; do j=1,3
      read(uf) bt(j,i)
    end do ;enddo
    b%tens = bt
  else
    write(*,'(A)') "Wrong size of box. Maybe that trr file is not supported by the present implementation."
    status = 1
    return
  end if
  if(np * 8 * 3 /= xsize) then
    write(*,'(A)') "Size of positions and nat does not match. Maybe that trr file is not supported by the present implementation."
    write(*,'("nat = ",I0," xsize = ", I0)') np, xsize 
    status = 1
    return
  end if
  b%volume = inverse_box(b%tens,boxstar)
  call calc_metric(b%tens, b%met)

  if(allocated(b%parts) .and. np > ubound(b%parts, 1)) deallocate(b%parts)
  if(.not. allocated(b%parts)) allocate(b%parts(np))
  b%n = np
  do i=1,b%n
    do j=1,3
      read(uf) pos(j)
    end do
    b%parts(i)%pos = matvec3_prod(boxstar,real(pos, 8))
  end do
  do i=1, vsize / 8
    read(uf) ld
  end do
  do i=1,fsize / 8
    read(uf) ld
  end do
end subroutine read_trr

end module trr

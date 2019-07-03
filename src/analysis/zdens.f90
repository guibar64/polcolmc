
module zdens
  ! Computes average particle density along the z-axis.
  use state, only: PcmcState
  implicit none

  real(8), allocatable :: histZdens(:,:,:)  !! Histogram of number of particle at a given plan XY. Array(z, fam, box)
  real(8) :: binZdens=1.D0                  !! space bin between two XY planes.
  integer :: nbinZdens=1000                 !! number of slices
  integer :: zdens_steps=0                  !! number of states sampled

  private :: histZdens, binZdens, nbinZdens, zdens_steps
  
  character(*), parameter, private :: suff_zdens="_zdens.dat"

contains

  subroutine zdens_init(st)
    !! Initialize calculations.
    type(PcmcState) :: st
    allocate(histZdens(nbinZdens, st%dist%nfam, st%ntotbox))
    histZdens = 0.0
    zdens_steps = 0
  end subroutine zdens_init

  subroutine zdens_update(st)
    !! Updates with current state.
    implicit none
    type(PcmcState) :: st
    integer :: i,k,nz
    real(8) :: Lz

    do k=1,st%ntotbox
      Lz = st%boxes(k)%tens(3,3)
      do i=1,st%boxes(k)%N
        nz = floor((st%boxes(k)%parts(i)%pos(3)*Lz)/binZdens) + 1
        if( nz <= nbinZDens ) histZdens(nz, st%boxes(k)%parts(i)%famille, k) = histZdens(nz, st%boxes(k)%parts(i)%famille, k) + 1.D0
      end do
    end do
    zdens_steps = zdens_steps + 1
  end subroutine zdens_update

  subroutine zdens_finish(st)
    !! Finalize calculations. `histZdens` is set to an average density.
    type(PcmcState) :: st
    histZdens = histZdens/(binZdens*(zdens_steps/st%inp%pZdens))
  end subroutine zdens_finish

  subroutine zdens_output(st, file_prefix)
    !! Prints densities to  files of the form `file_prefix`_`box_number`.dat
    type(PcmcState) :: st
    character(*) :: file_prefix
    integer :: k,i,j, ufo
    character(512) :: chardump
    do k=1,st%ntotbox
      write(chardump,'(A,"_box",I0,A)') file_prefix, k, suff_zdens 
      open(newunit=ufo, file=chardump)
      do i=1, nbinZdens
        write(ufo, '(ES16.7)', advance='NO') (real(i)-0.5D0)*binZdens
        do j=1,st%dist%nfam
          write(ufo, '(ES16.7)', advance='NO') histZdens(i, j,k)
        end do
        write(ufo,*)
      end do
    end do
    close(ufo)
  end subroutine zdens_output

  subroutine zdens_free()
    !! Frees arrays initialized in [[zdens_init]].
    deallocate(histZdens)
  end subroutine zdens_free
  
end module zdens

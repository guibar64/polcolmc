
module statcell
  !! Computes statistics on cell decomposition.
  use state, only: PcmcState
  implicit none

  integer :: n_statcell !! number of steps performed
  real(8),allocatable :: numofparts_statcell(:) !! average number of particles per cell
  real(8),allocatable :: var_numofparts(:)  !! variance of the number of particles per cell
  private :: n_statcell, numofparts_statcell, var_numofparts

  integer, allocatable, private :: numofparts_cells(:)
  
contains

  subroutine statcell_init(st)
    !! Initialize stats on cells.
    type(PcmcState) :: st
    allocate(numofparts_statcell(st%ntotbox),var_numofparts(st%ntotbox))

    n_statcell = 0
    numofparts_statcell(1:st%ntotbox) = 0.D0
    var_numofparts(1:st%ntotbox) = 0.D0
  end subroutine statcell_init

  subroutine statcell_update(st)
    !! Updates stats on cells.
    use celldec
    type(PcmcState) :: st
    integer :: k, ncellv
    
    n_statcell = n_statcell + 1
    do k=1,st%ntotbox
      call celldec_numofpart(st%clistdecomp(k),numofparts_cells, ncellv)
      numofparts_statcell(k) = numofparts_statcell(k) + st%boxes(k)%n/(1.D0*ncellv)
      var_numofparts(k) = var_numofparts(k) + sum(numofparts_cells(1:ncellv)*numofparts_cells(1:ncellv))&
           /ncellv
    end do
  end subroutine statcell_update

  subroutine statcell_finish(st)
    !! Finalize stat on cells.
    type(PcmcState) :: st
    if(n_statcell > 0) then
      numofparts_statcell(1:st%ntotbox) = numofparts_statcell(1:st%ntotbox)/n_statcell
      var_numofparts(1:st%ntotbox) = var_numofparts(1:st%ntotbox)/n_statcell - &
           numofparts_statcell(1:st%ntotbox)**2
    end if
  end subroutine statcell_finish

  subroutine statcell_log(st, log_unit)
    !! Prints stats on cells to a file.
    type(PcmcState) :: st
    integer, intent(in) :: log_unit !! file descriptor
    integer :: k
    if(n_statcell > 0) then
      do k=1,st%ntotbox
        write(log_unit,'(I8,G16.8,A,G16.8)') k, numofparts_statcell(k), "±", sqrt(var_numofparts(k))
      end do
    end if
  end subroutine statcell_log
  
  subroutine statcell_free()
    !! Free arrays used for computation. They are initializes by [[statcell_init]].
    deallocate(numofparts_statcell,var_numofparts)
  end subroutine statcell_free
  
end module statcell

module swapmap
  !! Computes a map of swaps between pairs of given families, ie the ratio
  !! between the number swaps done for a family pair and the total number of swaps.
  !! The same thing is also done for swaps between two simulation boxes. 
  !! Can be useful eg to check if swaps samples correctly.
  use state, only: PcmcState
  implicit none
  real(8),allocatable :: swap_map(:,:) !! Accepted swaps. Array(family, family)
  real(8),allocatable :: swap_map_box(:,:) !! Accepted swaps between simulation boxes.
  integer :: n_swap_map=0  !!  Number of sampled states.
  integer :: n_swap_map_box=0  !! Number of sampled states (box-box).

  ! Impurities
  integer,allocatable :: impur_stat_swap(:) !! Accepted swaps between "impurity" and a regular family. Array(family)
  integer :: n_impur_stat_swap=0 !!  Number of sampled states (inpurity).
  
contains

  subroutine swapmap_reset(st)
    !! Reset cumulative variables. Can be used before a new run.
    type(PcmcState) :: st
    n_swap_map = 0
    swap_map(1:st%dist%nfam,1:st%dist%nfam) = 0.d0
  end subroutine swapmap_reset
  
  subroutine swapmap_init(st)
    !! Initialize swapmap calculations.
    type(PcmcState) :: st
    allocate(swap_map(st%dist%nfam,st%dist%nfam))
    call swapmap_reset(st)
  end subroutine swapmap_init

  subroutine swapmap_update(i,j)
    !! Updates swapmap for a family pair (i,j).
    integer, intent(in) :: i, j
    n_swap_map = n_swap_map + 1
    swap_map(i, j) = swap_map(i, j) + 1.d0
  end subroutine swapmap_update

  subroutine swapmap_finish(st)
    !! Finalize swapmap calculations.
    type(PcmcState) :: st
    integer :: i,j
    if(n_swap_map > 0) then
      do i=1,st%dist%nfam
        do j=1,st%dist%nfam
          if(i<=j) swap_map(i,j)= 0.5D0*(swap_map(i,j) + swap_map(j,i))/(1.D0*n_swap_map)
        end do
      end do
    end if
  end subroutine swapmap_finish

  subroutine swapmap_output(st, file_name)
    !! Prints swapmap to a file.
    type(PcmcState) :: st
    character(*), intent(in) :: file_name
    integer :: i,j, ufo
    open(newunit=ufo,file=file_name)
    do i=1,st%dist%nfam-1
      do j=i+1,st%dist%nfam
        write(ufo,'(3ES16.7)') st%dist%rad(i), st%dist%rad(j), swap_map(i,j)
      end do
    end do
    close(ufo)
  end subroutine swapmap_output
  
  subroutine swapmap_free()
    !! Free swapmap array.
    deallocate(swap_map)
  end subroutine swapmap_free

  subroutine swapmap_box_reset(st)
    !! Reset box swapmap cumulative variables. Can be used before a new run.
    type(PcmcState) :: st
    n_swap_map_box = 0
    swap_map_box(1:st%dist%nfam,1:st%dist%nfam) = 0.d0
  end subroutine swapmap_box_reset
  
  subroutine swapmap_box_init(st)
    !! Initialize box swapmap calculations.
    type(PcmcState) :: st
    allocate(swap_map_box(st%dist%nfam,st%dist%nfam))
    call swapmap_box_reset(st)
  end subroutine swapmap_box_init

  subroutine swapmap_box_update(i,j)
    !! Updates box swapmap for a family pair (i,j).
    integer, intent(in) :: i, j
    n_swap_map_box = n_swap_map_box + 1
    swap_map_box(i, j) = swap_map_box(i, j) + 1.d0
  end subroutine swapmap_box_update

  subroutine swapmap_box_finish(st)
    !! Finalize box swapmap calculations.
    type(PcmcState) :: st
    integer :: i,j
    if(n_swap_map_box > 0) then
      do i=1,st%dist%nfam
        do j=1,st%dist%nfam
          if(i<=j) swap_map_box(i,j)= 0.5D0*(swap_map_box(i,j) + swap_map_box(j,i))/(1.D0*n_swap_map_box)
        end do
      end do
    end if
  end subroutine swapmap_box_finish

  subroutine swapmap_box_output(st, file_name)
    !! Prints box swapmap to a file.
    type(PcmcState) :: st
    character(*), intent(in) :: file_name
    integer :: i,j, ufo
    open(newunit=ufo,file=file_name)
    do i=1,st%dist%nfam-1
      do j=i+1,st%dist%nfam
        write(ufo,'(3ES16.7)') st%dist%rad(i),st%dist%rad(j),swap_map_box(i,j)
      end do
    end do
    close(ufo)
  end subroutine swapmap_box_output
  
  subroutine swapmap_box_free()
    !! Frees box swapmap array.
    deallocate(swap_map_box)
  end subroutine swapmap_box_free

  subroutine swapmap_impur_init(st)
    !! Initialize impur swapmap calculations.
    type(PcmcState) :: st
    allocate(impur_stat_swap(st%dist%nfam))
    n_impur_stat_swap = 0
    impur_stat_swap(1:st%dist%nfam) =0.D0
  end subroutine swapmap_impur_init

  subroutine swapmap_impur_update(fam)
    !! Update impur swapmap for family `fam`.
    integer, intent(in) :: fam
    n_impur_stat_swap=n_impur_stat_swap+1
    impur_stat_swap(fam) = impur_stat_swap(fam) + 1
  end subroutine swapmap_impur_update

  subroutine swapmap_impur_log(st, log_unit)
    !! Prints information about impur swaps.
    type(PcmcState) :: st
    integer, intent(in) :: log_unit
    integer :: i
    write(log_unit,'(A)') "Impur swaps : families stats"
    do i=1,st%dist%nfam
      write(log_unit,'(2G16.8)') st%dist%rad(i), (1.d0*impur_stat_swap(i))/(1.d0*n_impur_stat_swap)
    end do
  end subroutine swapmap_impur_log
  
  subroutine swapmap_impur_free()
    !! Frees impur swapmap array.
    deallocate(impur_stat_swap)
  end subroutine swapmap_impur_free
  
end module swapmap

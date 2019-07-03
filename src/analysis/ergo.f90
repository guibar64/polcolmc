module ergo
  !! Ergodicity analysis based on particle energy correlations
  use state
  implicit none
  
  type NrgErgo
    !! for energy ergo analysis
    real(8), allocatable :: cumepp(:) !! cumulative step particle energies
    integer :: nsteps
    real(8), allocatable :: mean(:), sigma(:)  !! mean, variance per particle family
  end type NrgErgo
contains

  subroutine nrgergo_init(nea, st)
    !! Initialize analysis.
    type(NrgErgo) :: nea(*)
    type(PcmcState), intent(in) :: st
    integer :: k
    do k=1,st%ntotbox
      if(allocated(nea(k)%cumepp)) deallocate(nea(k)%cumepp, nea(k)%mean, nea(k)%sigma)
       allocate(nea(k)%cumepp(ubound(st%boxes(k)%parts, 1)), nea(k)%mean(st%dist%nfam),&
        nea(k)%sigma(st%dist%nfam))
      nea(k)%cumepp = 0.d0
      nea(k)%nsteps = 0
      nea(k)%mean = 0.d0
      nea(k)%sigma = 0.d0
    end do
  end subroutine nrgergo_init

  subroutine nrgergo_update(nea, b)
    !! Updates `ǹea` with a new state
    implicit none
    type(NrgErgo) :: nea
    type(Box), intent(in) :: b
    integer :: i, fam
    nea%nsteps = nea%nsteps + 1
    nea%mean = 0.d0
    nea%sigma = 0.d0
    do i=1,b%n
      fam = b%parts(i)%famille
      nea%mean(fam) = nea%mean(fam) + nea%cumepp(i)
      nea%sigma(fam) = nea%sigma(fam) + nea%cumepp(i)**2
    end do
    do fam=1,ubound(b%fam,1)
      if(b%fam(fam) > 0) then
        nea%mean(fam) = nea%mean(fam) / (b%fam(fam)*nea%nsteps)
        nea%sigma(fam) = nea%sigma(fam) / (b%fam(fam)*real(nea%nsteps,8)**2) - nea%mean(fam)**2
      end if
    end do
  end subroutine nrgergo_update

  subroutine nrgergo_write(nea, step, unit)
    !! Prints current averages to `unit`. 
    implicit none
    type(NrgErgo), intent(in) :: nea
    integer, intent(in) :: unit !! output file unit
    integer :: step !! current step (should begin at 0) 
    integer :: j
    write(unit,'(I9)', advance='no') step
    do j=1, ubound(nea%sigma, 1)
      write(unit,'(g16.7)', advance='no') nea%sigma(j)  
    end do
    do j=1, ubound(nea%mean, 1)
      write(unit,'(g16.7)', advance='no') nea%mean(j)
    end do
    write(unit, *)    
  end subroutine nrgergo_write

  subroutine nrgergo_free(nea)
    type(NrgErgo) :: nea
    deallocate(nea%cumepp, nea%mean, nea%sigma)
  end subroutine nrgergo_free

end module ergo
module chempot
  !! Calculations of chemical potentials
  use state, only: PcmcState, Particule, NrgRec, kT
  implicit none

  ! ChemPot calculations globals
  logical :: chempwid_initialized = .false.   !! toggles the computation or not of the chemical otentials
  ! Widom-like way
  integer,allocatable :: chempwid_nsample(:,:) !! number of samples (widom). Array(box, family]
  real(8),allocatable :: chempwid_bmu(:,:)     !! excess chemical potential (cumulated during simulation). Array(box, family)
  real(8),allocatable :: chempwid_bmu_tot(:,:) !! chemical potential. Array(box, family)
  ! With "Rosenbluth" moves
  integer,allocatable :: chempros_nsample(:,:)  !! number of samples (Rosenbluth). Array(box, family)
  real(8),allocatable :: chempros_bmu(:,:)      !! excess chemical potential (cumulated during simulation). Array(box, family)
  real(8),allocatable :: chempros_bmu_tot(:,:)  !! chemical potential. Array(box, family)

contains

  subroutine chempwid_reset(st)
    !! Reset all cumulative quantitive. Useful eg when doing a new run.
    type(PcmcState) :: st
    chempwid_nsample(1:st%dist%nfam,1:st%ntotbox) = 0
    chempwid_bmu(1:st%dist%nfam,1:st%ntotbox) = 0.d0
    chempros_nsample(1:st%dist%nfam,1:st%ntotbox) = 0
    chempros_bmu(1:st%dist%nfam,1:st%ntotbox) = 0.d0
  end subroutine chempwid_reset
  
  subroutine chempwid_init(st)
    !! Initialize chempot calculationss
    type(PcmcState) :: st
    allocate(chempwid_nsample(st%dist%nfam,st%ntotbox),chempros_nsample(st%dist%nfam,st%ntotbox))
    allocate(chempwid_bmu(st%dist%nfam,st%ntotbox),chempros_bmu(st%dist%nfam,st%ntotbox))
    allocate(chempwid_bmu_tot(st%dist%nfam,st%ntotbox),chempros_bmu_tot(st%dist%nfam,st%ntotbox))
    chempwid_initialized = .true.
    call chempwid_reset(st)
  end subroutine chempwid_init

  subroutine chempwid_update(st, nr, nsamples)
    !! sample excess μ with Widom insertions
    use moves, only: add_part_box, remove_part_box
    use random2
    use nrg
    type(PcmcState), intent(in) :: st
    type(NrgRoutines), intent(in) :: nr
    type(Particule) :: part
    integer, intent(in) :: nsamples !! number of samples
    integer :: i,k, ip, s, j, d
    type(NrgRec) :: e
    do k=1,st%ntotbox
      do i=1,st%boxes(k)%n
        j = st%boxes(k)%parts(i)%famille
        ip = st%boxes(k)%n + 1
        ! Transfer properties
        part = st%boxes(k)%parts(i)
        do s=1,nsamples
          do d=1,3
            part%pos(d) = ran2(st%rng_an)
          end do
          call add_part_box(st, k, part)
          e = nr%energy_1p(st, st%boxes(k), st%boxes(k)%parts(ip))
          if(e%noverlap) chempwid_bmu(j,k) = chempwid_bmu(j,k) + exp(-e%ep/kT)
          call remove_part_box(st, k, ip)
        end do
        chempwid_nsample(j, k) = chempwid_nsample(j, k) + nsamples
      end do
    end do

  end subroutine

  subroutine chempwid_finish(st)
    !! To be called at the end of a run. chemp*_bmu(i,k) is reduced and chemp*_bmu_tot(i,k) is set.
    type(PcmcState) :: st
    integer :: k,i
    do k=1,st%ntotbox
      do i=1,st%dist%nfam
        if(chempwid_nsample(i,k)>0 .and. st%res(k)%popm(i)>=1.D0 .and. chempwid_bmu(i,k)>0.D0) then
          chempwid_bmu(i,k) = - log(chempwid_bmu(i,k)/chempwid_nsample(i,k))
          chempwid_bmu_tot(i,k) = chempwid_bmu(i,k) + log(st%res(k)%popm(i)/st%res(k)%volm)
        else
          chempwid_bmu(i,k) = 0.D0
          chempwid_bmu_tot(i,k) = 0.D0
        end if
        if(chempros_nsample(i,k)>0 .and. st%res(k)%popm(i)>=1.D0 .and. chempros_bmu(i,k)>0.D0) then
          chempros_bmu(i,k) = - log(chempros_bmu(i,k)/chempros_nsample(i,k))
          chempros_bmu_tot(i,k) = chempros_bmu(i,k) + log(st%res(k)%popm(i)/st%res(k)%volm)
        else
          chempros_bmu(i,k) = 0.D0
          chempros_bmu_tot(i,k) = 0.D0
        end if
      end do
    end do
  end subroutine chempwid_finish

  subroutine chempwid_output(st, file_name)
    !! Prints the chem*_bmu_[_tot] to a file.
    type(PcmcState) :: st
    character(*), intent(in) :: file_name
    integer :: j,k, ufo
    open(newunit=ufo,file= file_name)
     do j=1,st%dist%nfam
       write(ufo,'(ES16.7)',advance='no') st%dist%rad(j)
       do k=1,st%ntotbox
         write(ufo,'(4ES16.7)',advance='no') chempwid_bmu_tot(j,k),chempwid_bmu(j,k),&
              chempros_bmu_tot(j,k),chempros_bmu(j,k)
       end do
       write(ufo,*)
     end do
     close(ufo)
  end subroutine chempwid_output
  
  subroutine chempwid_free()
    !! Free all arrays. Should be called when they are no longer needed.
    deallocate(chempwid_nsample, chempros_nsample)
    deallocate(chempwid_bmu, chempros_bmu)
    deallocate(chempwid_bmu_tot, chempros_bmu_tot)
    chempwid_initialized = .false.
  end subroutine chempwid_free

end module chempot

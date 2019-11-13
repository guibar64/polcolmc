module misc
!! Miscellaneous utilities
implicit none


  ! need to make dist2_* private to avoid re-declarations
  include "dist2_privs.f90"

contains

  ! dist2_* routines
  include "dist2_mi.f90"
  
  real(8) function dist2_mi(pos,posout, metricG,m_box)
    !! Returns qquare of the separation between the minimum image of two particles
    !! `posout` is set to the vector between the images.
    implicit none
    real(8),intent(in) :: pos(3) !! "real" separation vector
    real(8),intent(out) :: posout(3) !! separation vector
    real(8),intent(in) :: metricG(3,3) !! box metric
    integer, intent(in) :: m_box !! box type
    select case(m_box)
    case(0)
      dist2_mi = dist2_mi_generic(pos,posout, metricG)
    case(1)
      dist2_mi = dist2_mi_cubic(pos,posout, metricG)
    case(2)
      dist2_mi = dist2_mi_cubicXY(pos,posout, metricG)
    case default
      dist2_mi = dist2_mi_generic(pos,posout, metricG)
    end select
  end function dist2_mi

  subroutine search_neighbors(st,rad,part,ibox,neighn,neighlist,neighmax)
    !! Does a neigbor search around particle part on box `ibox`.
    use state, only: PcmcState, Box
    use celldec, only: Particule
    type(PcmcState) :: st
    type(Particule), intent(in) :: part
    integer,intent(in) :: neighmax,ibox
    integer, intent(out) :: neighlist(neighmax),neighn
    real(8), intent(in) :: rad

    integer :: i
    real(8) :: r2,pos(3), prepos(3)

    associate(labox => st%boxes(ibox))

    ! TODO: optimize with cell decomposition when fit
    neighn = 0
    do i=1, labox%n
      if(part%numero /= labox%parts(i)%numero) then
        prepos = part%pos - labox%parts(i)%pos
        r2 = dist2_mi(prepos,pos, labox%Met,st%mode_box)
        if(r2<rad) then
          neighn=neighn+1
          neighlist(neighn)=labox%parts(i)%numero
          if(neighn>=neighmax) exit
        end if
      end if
    end do
    end associate
  end subroutine search_neighbors

  logical function random_placement(n,Partics,graine,metricG,m_box)
    !! Place a set of particles in a box at random positions
    !! so that the particles do not overlap.
    use celldec, only: Particule
    use random2
    implicit none
    integer,intent(in) :: n !! number of particles
    integer,intent(in) :: m_box !! box type
    integer,intent(inout) :: graine !! random seed
    type(Particule),intent(inout) :: Partics(n) !! particles
    real(8) :: metricG(3,3) !! box metric

    logical :: istok
    integer,parameter :: maxiter=10000
    integer :: i,j,it
    real(8) :: dump
    real(8) :: posrel(3), prepos(3)
    type(Ran2State) :: randg
    call ran2_put_seed(randg, graine)

    ! Assuming that particles are ordered according to their decreasing size

    do i=n,1,-1
      do it=1,maxiter
        istok=.true.
        do j=1,3
          Partics(i)%pos(j)=RAN2(randg)
        enddo
        do j=i+1,n,1
          prepos=Partics(i)%pos-Partics(j)%pos
          dump=dist2_mi(prepos, posrel,metricG,m_box)
          if(dump<=(Partics(i)%rayon+Partics(j)%rayon)**2) then
            istok=.false.
            exit
          end if
        end do
        if(istok) exit
      end do
      if(it>maxiter) then
        random_placement=.false.
        return
      endif
    end do

    random_placement=.true.

  end function random_placement

  logical function isnot_overlapped(n,parts,metricG, ifault, jfault,m_box) result(over)
    !! Checks if there is overlaps in a particle state.
    !! On overlap, the first overlapping pair found is returned to (`ifault`, `jfault`)
    use celldec, only: Particule
    implicit none
    integer,intent(in) :: n !! number of particles
    integer,intent(in) :: m_box !! box type
    integer,intent(out) :: ifault, jfault
    real(8),intent(in) :: metricG(3,3) !! box metric
    type(Particule),intent(in),dimension(n) :: parts !! particles

    integer :: j,i
    real(8) :: r2
    real(8) :: prepos(3), relpos(3)

    ! Could benefit from cell decomposition. Though it is currently only used if something went wrong.
    do i=1,n
      do j=i+1,n
        if(i == j) cycle
        prepos=parts(i)%pos-parts(j)%pos
        r2=dist2_mi(prepos,relpos,metricG,m_box)
        if(r2<(parts(i)%rayon+parts(j)%rayon)**2) then
          ifault = i; jfault = j
          over=.false.
          return
        endif
      end do
    end do

    over=.true.
  end function isnot_overlapped

  subroutine short_initialisation(st,fichdist,fichpar)
    !! Only initialize distribution and input parameters
    use readparam
    use state
    implicit none
    type(PcmcState) :: st
    character(*), intent(in) :: fichdist !! distribution file
    character(*), intent(in) :: fichpar  !! configuration file
    integer :: stat

    ! Read distribution file
    call init_distribution(st%dist,fichdist, stat)
    if(stat<0) then
      write(*,'(A)') "Error: distribution '",trim(fichdist),"' cannot be found."
      stop 1
    endif
    st%Nbuff = sum(st%dist%pop)

    ! Read configuration file
    call init_param_list(st%config)
    call read_cfg_file(st,fichpar)
    if (stat /= 0) then
      write(*,'(A,A,A)') "Error: Configuration file '", trim(fichpar),"' cannot be found."
      stop 1
    end if
    call configuration_finish(st)
  end subroutine short_initialisation

  subroutine  pos_max1(st, potent,i,j,rmin,rmax,max)
    !! Finds the first maximum of a potential `potent` between families `i` and `j`,
    !! from separation `rmin`.
    !! `rmax` is set to the separation of the maximum and `max` to the value of
    !! the max potential  
    use state, only: PcmcState
    implicit none
    external :: potent
    type(PcmcState) :: st
    real(8) :: potent
    integer,intent(in) :: i,j
    real(8), intent(in) :: rmin
    real(8),intent(out) :: rmax,max

    real(8) :: r0,r1,r2,dr,pot0,pot1,pot2,tol
    integer :: ii,n_it
    parameter(tol=1.e-3_8,dr=0.1_8,n_it=200)
    r0=rmin
    pot0=potent(r0,i,j,st)

    ! Finds a maximum by evaluating at constant interval dr. This is obviously wrong in general...
    do ii=1,n_it
      r1=r0+dr
      pot1=potent(r1,i,j,st)
      if(pot1<pot0) exit

      r0=r1
      pot0=pot1
    enddo
    ! Bisecting now.
    do ii=1,n_it
      r2=r0+(r1-r0)/2
      pot2=potent(r2,i,j,st)
      if((pot2-pot0)<(pot2-pot1)) then
        r1=r2
        pot1=pot2
      else
        r0=r2
        pot0=pot2
      endif
      if(abs((pot1-pot0))<tol) exit
    enddo
    rmax=r0+(r1-r0)/2
    max=potent(r0+(r1-r0)/2,i,j,st)

  end subroutine pos_max1

end module misc



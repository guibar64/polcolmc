module correls
  !! Calculates radial distribution functions (RDF) and effective structure factors (S(q)).
  use celldec, only: Particule
  implicit none

  type CorrelsRec
     !! RDF calculation type
     integer :: ng  !! number of separations
     integer ::  nr !! number of families
     integer :: step_calc !! number of steps done
     real(8) :: dr !! space bin
     real(8) :: rmin !! minimum radius
     real(8),allocatable :: g(:,:) !! RDF. Array(ng, nr*(nr+1)/2)
     real(8),allocatable :: gmoy(:) !! mean (or total) RDF. Array(ng)
     integer(8),allocatable :: nder(:,:) !! cumulant for computation
  end type CorrelsRec

  type StructqRec
    !! structure factors type. There are two versions: normal and "full"
    !! Full uses the q-values of the DFT of the RDFs.
    !! The normal version uses input [`qmin`, `q_max`] and bin `dq`. 
    integer :: ns !! number of q-values
    integer :: nsfull !! number of q-values (full)
    real(8) :: dq !! reciprocal space bin
    real(8) :: qmin !! minimum q value.
    real(8),allocatable :: S(:) !! S(q)
    real(8),allocatable :: q(:) !! q-values 
    real(8),allocatable :: qfull(:) !! q-values (full)
    real(8),allocatable :: Sfull(:) !! S(q) (full)
    real(8),allocatable :: formf(:) !! form factor (full)
  end type StructqRec


  real(8),parameter :: PI=3.1415926535897931D0

  ! interface
  !   subroutine correls_new_step_type(corrs,metricG,n,centres)
  !     !! Prototype of a function that update the RDF histograms.
  !     import
  !     type(CorrelsRec),intent(inout) :: corrs !! RDF calculator
  !     integer,intent(in) :: n !! number of particules
  !     real(8),intent(in),dimension(3,3) :: MetricG !! box metric
  !     type(Particule),intent(in) :: Centres(*) !! particle array
  !   end subroutine correls_new_step_type
  ! end interface
  ! procedure(correls_new_step_type),pointer :: correls_new_step !! Do a new step 
  
  private :: PI,Indic_

  private :: correls_new_step_def, correls_new_step_cubic
  
  include "dist2_privs.f90"

contains

! Not DRY.
include "dist2_mi.f90"

subroutine correls_init(corrs,n_r,d_r,r_min,r_max)
  !! Initialize an RDF calculation.
  implicit none
  type(CorrelsRec),intent(out) :: corrs
  integer,intent(in) :: n_r !! number of families
  real(8),intent(in) :: d_r !! space bin
  real(8),intent(in) :: r_min !! minimum radius
  real(8),intent(in) :: r_max !! maxmimum radius

  integer :: ng

  corrs%dr=d_r
  ng=1+floor((r_max-r_min)/d_r)
  corrs%ng=ng
  corrs%nr=n_r
  corrs%rmin=r_min

  allocate(corrs%nder(ng,(n_r*(n_r+1))/2))

  corrs%step_calc=0
  corrs%nder=0

end subroutine correls_init

subroutine correls_reset(corrs)
  !! Reset RDF calculations. Can be called with a new run.
  implicit none
  type(CorrelsRec),intent(inout) :: corrs
  corrs%step_calc=0
  corrs%nder=0
  
end subroutine correls_reset

subroutine correls_new_step(corrs,metricG,n,centres,mode)
  !! Do a new RDF step
  type(CorrelsRec),intent(inout) :: corrs !! RDF calculator
  integer,intent(in) :: n !! number of particules
  real(8),intent(in),dimension(3,3) :: MetricG !! box metric
  type(Particule),intent(in) :: Centres(*) !! particle array
  integer, intent(in) :: mode              !! box mode

  select case(mode)
  case(0) ! general periodic box
     call correls_new_step_def(corrs,metricG,n,centres)
  case(1) ! cubic periodic box
     call correls_new_step_cubic(corrs,metricG,n,centres)
  case(2) ! Z Hard Walls
     ! cheating, as non-PBC are not compatible with homogeneous g(r)'s
     call correls_new_step_cubic(corrs,metricG,n,centres)  
  case default
     call correls_new_step_def(corrs,metricG,n,centres)
  end select
end subroutine correls_new_step

subroutine correls_end_step(corrs,vol,pop)
  !! Compute RDF after a run. Sets the fields `g` of `corrs`.
  implicit none
  type(CorrelsRec),intent(inout) :: corrs
  real(8),intent(in) :: pop(corrs%nr) !! family populations
  real(8),intent(in) :: vol   !! volume of the simulation box
  integer :: i,j,k,ind
  real(8) :: dump

  real(8),parameter :: EPSI=1.1D0
  integer :: nr,ng,step_calc
  real(8) :: dr

  nr=corrs%nr
  ng=corrs%ng
  dr=corrs%dr
  step_calc=corrs%step_calc

  allocate(corrs%g(ng,(nr*(nr+1))/2))
  if(step_calc>0) then

     dump=vol*((1.D0/dr)**3)/(4*PI)
     corrs%g = (dump*corrs%nder)/step_calc

     do i=1,nr
        if(pop(i)>EPSI) then
           do k=1,ng
              corrs%g(k,i)=corrs%g(k,i)/(k*(k-1)+0.333333333D0)/(pop(i)*(pop(i)-1.d0))
           end do
        else
           corrs%g(1:ng,i)=0.D0
        end if
     end do
     do i=1,nr
        do j=i+1,nr
           ind=Indic_(i,j,nr)
           if(pop(i)>EPSI.and.pop(j)>EPSI) then
              do k=1,ng
                 corrs%g(k,ind)=0.5D0*corrs%g(k,ind)/(k*(k-1)+0.333333333D0)/(pop(i)*pop(j))
              enddo
           else
              corrs%g(1:ng,ind)=0.D0
           end if
        end do
     end do
  else
     corrs%g(:,:)=0.D0
  end if

end subroutine correls_end_step

subroutine correls_free(corrs)
  !! Free arrays of `corrs` (allocated by [[correls_init]])
  implicit none
  type(CorrelsRec),intent(inout) :: corrs
  if(allocated(corrs%nder)) deallocate(corrs%nder)
  if(allocated(corrs%g)) deallocate(corrs%g)

end subroutine correls_free


subroutine correls_reduit_g(corrs,poprad)
  !! Compute mean RDF. Sets the fields `gmoy` of `corrs`.
  implicit none
  type(CorrelsRec),intent(inout) :: corrs
  real(8),intent(in) :: poprad(corrs%nr) !! family populations

  real(8),parameter :: EPSI=1.D-6
  integer :: i,j,nr
  real(8) :: totpop
  
  nr=corrs%nr

  Totpop=sum(poprad(1:nr))
  
  allocate(corrs%gmoy(corrs%ng))
  corrs%gmoy(:)=0.D0
  if(totpop>=EPSI) then
     do j=1,nr
        if(Poprad(j)>EPSI) then
           corrs%gmoy(:)=corrs%gmoy(:)+corrs%g(:,j)*(poprad(j)*(poprad(j)-1.d0))
        end if
     enddo

     do i=1,nr ; do j=i+1,nr
        if(poprad(i)>EPSI.and.poprad(j)>EPSI) then
           corrs%gmoy(:)=corrs%gmoy(:)+2*corrs%g(:,((2*nr-i-1)*i)/2+j)*(poprad(i)*poprad(j))
        end if
     enddo;enddo
  end if
  corrs%gmoy(:)=corrs%gmoy(:)/(totPop*(totpop-1.d0))

end subroutine correls_reduit_g


integer function indic_(i,j,n)
   !! Family indices to index in a flat sequence of family-family pairs 
   !! (1,1),..,(n,n),(1,2),(1,3),..,(2,3),..(n-1,n)
   integer,intent(in) :: i,j
   integer,intent(in) :: n !! number of families
   if(i==j) then
      indic_=i
   else
      if(i<j) then
         indic_=((2*n-i-1)*i)/2+j
      else
         indic_=((2*n-j-1)*j)/2+i
      endif
   endif
 end function indic_

! Preprocessor hack to generate specialized routines (according to box type)

! Generic Box
#define DIST2_MI dist2_mi_generic
#include "correls_default.f90"
#undef DIST2_MI

! Cubic Box
#define DIST2_MI dist2_mi_cubic
#define CORRELS_NEW_STEP_DEF correls_new_step_cubic
#include "correls_default.f90"
#undef DIST2_MI
#undef CORRELS_NEW_STEP_DEF


subroutine structq_free(sq)
  !! Frees arrays of `sq`. They are initialized by [[SofqfromGder]].
  type(StructqRec) :: sq
  deallocate(sq%q,sq%S,sq%qfull,sq%Sfull,sq%formf)    
end subroutine

subroutine SofqfromGder(alls,corrs,volume,Poprad,Radius,qmin,qmax,dq)
  !! Calculates S(q) from the RDFs. Results to `alls`.
  use struct_fact
  implicit none
  type(CorrelsRec),intent(in) :: corrs !! RDFs
  type(StructqRec), intent(inout) :: alls  !! S(q) results
  real(8),intent(in) :: volume !! volume of the simulation box
  real(8),intent(in) :: Radius(corrs%nr) !! Family radii
  real(8),intent(in) :: Poprad(corrs%nr) !! family populations 

  integer :: i,nr,ng,ns,nfft
  real(8) :: dq,qmin,qmax,dr,totpop,rho, phi
  real(8),allocatable, target :: normS(:),ss(:), qfft(:)
  
  nr=corrs%nr
  ng=corrs%ng
  dr=corrs%dr

  ns=floor((qmax-qmin)/dq)
  alls%ns=ns
  allocate(alls%S(ns), alls%q(ns))
  do i=1,ns
     alls%q(i) = qmin+ (i-1)*dq
  end do
  
  nfft=ng
  allocate(normS(nfft),ss(nfft),qfft(nfft))

  totpop = sum(poprad(1:nr))
  rho = totpop/volume
 
  phi = sum(poprad*4./3.*PI*Radius**3)/volume

  if(corrs%step_calc > 0) then
    call calc_sdeq(dr, ng, qfft , ss , normS, ns, alls%q, alls%s, (nr*(nr+1))/2, corrs%g, nr, poprad, radius, phi)
  end if
  call move_alloc(qfft, alls%qfull)
  call move_alloc(ss, alls%Sfull)
  call move_alloc(norms, alls%formf)
  alls%nsfull = ng
  
end subroutine

end module CORRELS


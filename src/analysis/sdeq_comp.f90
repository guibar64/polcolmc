module sdeq_comp
!! Calculates directly the struture factors as function of q-vectors.
use celldec, only: Particule
implicit none

type Qvector
  !! Wave vector
  real(8) :: comp(3) !! coordinates
  real(8) :: norm    !! norm
end type Qvector

real(8),parameter :: PI=3.1415926535897931D0
integer :: fusc
integer :: nq    !! number of q vectors
integer :: nstep !! number of computed steps
type(Qvector) :: q  !! q vectors
real(8) :: bigs !! S(q)
real(8) :: norm  !! form factor
allocatable :: bigs(:),q(:),norm(:)

private :: q,nq,nstep,formfact,PI,fusc,norm,Qvector

Contains

subroutine sdeq_comp_init(fichq)
  !! Initialize calculations
  character(*),intent(in) :: fichq

  logical :: fichestla
  integer :: i

  inquire(file=fichq,exist=fichestla)
  if(.not.fichestla) then
    write(*,'(A)') "Cannot open wave number file."
    nq=0
    stop 2
  else
    open(newunit=fusc,file=fichq,status='old')
    read(fusc,*) nq
    allocate(bigs(nq),q(nq),norm(nq))
    do i=1,nq
      read(fusc,*) q(i)%comp
      q(i)%norm=sqrt(dot_product(q(i)%comp,q(i)%comp))
    end do
  end if

  close(fusc)

  bigs=0.D0
  nstep=0

end subroutine sdeq_comp_init

subroutine sdeq_comp_reset()
  !! reset calculations. Can be used for a new run.
  bigs=0.d0
  nstep=0

end subroutine sdeq_comp_reset

subroutine sdeq_comp_new_step(parts,n)
  !! Updates with current state.
  integer,intent(in) :: n !! number of particles
  type(Particule),intent(in) :: parts(n) !! particles

  integer :: i,j
  real(8) :: x
  complex(8) :: phi_s(nq)

  phi_s(:)=(0.d0,0.d0)
  do i=1,n
    do j=1,nq
      x=dot_product(parts(i)%pos,q(j)%comp)
      phi_s(j)=phi_s(j)+(formfact(q(j)%norm*parts(i)%rayon))*exp(cmplx(0.d0,x,kind=8))
    end do
  end do

  bigs(1:nq)=bigs(1:nq)+real(phi_s*conjg(phi_s), 8)
  nstep=nstep+1

end subroutine sdeq_comp_new_step

subroutine sdeq_comp_output(fichout)
  !! Prints results (S(q), P(q)) to a file.
  implicit none
  character(*),intent(in) :: fichout !! path to output file.

  integer :: j, ufo

  open(newunit=ufo,file=fichout)

  do j=1,nq
    write(fusc,'(5g16.7)') q(j)%comp,bigs(j)/norm(j),norm(j)
  end do
  close(ufo)
end subroutine sdeq_comp_output

subroutine sdeq_comp_end_step(Parts,n)
  !! Finalize calculations of S(q).
  implicit none
  integer,intent(in) :: n
  type(Particule),intent(in) :: Parts(n)

  integer :: i,j
  real(8) :: vol2

  bigs=bigs/(nstep)

  norm=0.D0
  vol2=0.D0
  do i=1,n
    do j=1,nq
      norm(j)=norm(j)+formfact(q(j)%norm*Parts(i)%rayon)**2
    end do
    vol2=vol2+(4.d0/3.d0*PI*(Parts(i)%rayon)**3)**2
  end do
  
  do j=1,nq
    bigs(j)=bigs(j)/(vol2*(q(j)%norm)**6)
    norm(j)=norm(j)/(vol2*(q(j)%norm)**6)
  end do

end subroutine sdeq_comp_end_step

real(8) function formfact(x)
  !! Form factor of an homogeneous sphere.
  real(8) :: x

  formfact=3*(sin(x)-x*cos(x))
end function formfact

end module sdeq_comp

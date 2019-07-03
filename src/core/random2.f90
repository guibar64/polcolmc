module random2
  !! A pseudo-random number generator.
  implicit none
  integer,private ::IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  real(8),private :: AM,EPS,RNMX
  parameter (IM1=2147483563,IM2=2147483399,AM=1.0/IM1,IMM1=IM1-1,&
  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
  IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-9,RNMX=1.-EPS)

 type Ran2State
    !! RNG state
    integer :: idum
    integer :: idum2,iv(NTAB),iy
 end type Ran2State
 
contains

subroutine ran2_reset(s)
  !! Reset ran2 RNG
  type(Ran2State) :: s
  s%idum2=123456789
  s%iv = 0
  s%iy = 0
end subroutine ran2_reset

subroutine ran2_put_seed(s, idum)
  !! Seeds ran2 RNG 
  integer, intent(in) :: idum
  type(Ran2State) :: s
  integer :: j,k
  s%idum=max(abs(idum),1)
  s%idum2=s%idum
  do j=NTAB+8,1,-1
     k=s%idum/IQ1
     s%idum=IA1*(s%idum-k*IQ1)-k*IR1
     if (s%idum.lt.0) s%idum=s%idum+IM1
     if (j.le.NTAB) s%iv(j) = s%idum
  enddo
  s%iy = s%iv(1)
end subroutine ran2_put_seed

  
real(8) function ran2(s)
  !! Returns a random real in [0,1[.
  !! Repeated calls form a random sequence. 
  type(Ran2State) :: s

  integer :: j,k
  k=s%idum/IQ1
  s%idum=IA1*(s%idum-k*IQ1)-k*IR1
  if (s%idum.lt.0) s%idum = s%idum+IM1
  k=s%idum2/IQ2
  s%idum2=IA2*(s%idum2-k*IQ2)-k*IR2
  if (s%idum2.lt.0) s%idum2=s%idum2+IM2
  j=1+s%iy/NDIV
  s%iy=s%iv(j)-s%idum2
  s%iv(j)=s%idum
  if(s%iy.lt.1) s%iy=s%iy+IMM1
  ran2=min(AM*s%iy,RNMX)
end function ran2

logical function metropolis(dx,rng)
  !! Accepts or rejects by the metropolis rule.
  !! Returns `.true.` if min(1, exp(-dx)) >= (random number),
  !! .false. otherwise.
  type(Ran2State) :: rng   !! RNG
  real(8),intent(in) :: dx !! can be seen as ΔE/kT
  real(8) :: dump
  if(dx<=0) then
     metropolis=.true.
     return
  endif
  dump=exp(-dx)
  if(ran2(rng)<dump) then 
     metropolis=.true.
  else
     metropolis=.false.
  endif
end function metropolis

end module random2

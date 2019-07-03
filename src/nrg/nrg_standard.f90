function NRG_1P(st,b,par)
  implicit none
  type(NrgRec) :: NRG_1P
  type(PcmcState), intent(in) :: st
  type(Box), intent(in),target :: b
  type(Particule),intent(in),target :: par

  integer :: i,n_r
  real(8) :: e,r2,phs1,phs2,pvir
  real(8) :: prepos(3), relpos(3),rc2
  type(Particule), pointer :: curr_par
  n_r = st%dist%nfam
  rc2 = st%rcrc
  e=0.D0
  phs1=0.D0
  pvir=0.D0
  phs2=0.D0

  !$omp parallel do default(shared),private(i,r2,curr_par,prepos,relpos)&
  !$omp reduction(+:e)
  do i=1,b%n

     curr_par => b%parts(i)
     
     if(.not.associated(curr_par,par)) then
        prepos=par%pos-curr_par%pos
        r2=DIST2_MI(prepos, relpos,b%met)

        if(r2<(par%rayon+curr_par%rayon)**2) then
          NRG_1P%noverlap=.false.
#ifndef _OPENMP
          return
#endif
        elseif(r2<rc2) then
           e = e + curr_par%ech*potent_r2(r2, st%ndat)
        end if
     end if
  end do
  NRG_1P%noverlap = .true.
  NRG_1P%ep = st%ndat%elpot_pref * par%ech * e
  
end function NRG_1P

function NRG_BOX(st,b,press)
  implicit none
  type(NrgRec) :: NRG_BOX
  type(PcmcState), intent(in) :: st
  type(Box), intent(in),target :: b
  type(Pressure),intent(out),optional :: press

  logical :: pressons
  integer :: j,i, n_r
  real(8) :: dump,r,e,ei,r2,phs1,phs2,pvir,pviri
  real(8) :: prepos(3), relpos(3),rc2, dr_press
  type(Particule), pointer :: curr_par,par
  pressons = present(press)
  dr_press = st%inp%dr_press
  n_r = st%dist%nfam
  rc2 = st%rcrc
  e=0.D0
  phs1=0.D0
  pvir=0.D0
  phs2=0.D0

  !$omp parallel do default(shared),private(i,j,r2,r,par,curr_par,prepos,relpos,dump,ei,pviri)& 
  !$omp reduction(+:e,phs2,phs1,pvir)
  do i=1,b%n
     ei = 0.d0
     pviri = 0.d0
     par => b%Parts(i)
     do j=i+1,b%n
        curr_par => b%Parts(j)
     
        if(.not.associated(curr_par,par)) then
           prepos=par%pos-curr_par%pos
           r2=DIST2_MI(prepos,relpos,b%met)
           if(r2<(par%rayon+curr_par%rayon)**2) then
              NRG_BOX%noverlap=.false.
#ifndef _OPENMP
              return
#endif
           elseif(r2<rc2) then
              if(pressons) then
                 r=sqrt(R2)
                 ei = ei + curr_par%ech*potent_r(r, st%ndat)
                 pviri = pviri + curr_par%ech*r*force_r(r, st%ndat)
              
                 dump=par%rayon+curr_par%rayon
                 if(r<dump+dr_press) then
                    phs2=phs2+dump/(2*dr_press)
                    phs1=phs1+dump/dr_press
                 elseif(r<dump+2*dr_press) then
                    phs2=phs2+dump/dr_press
                 endif
              else
                 ei = ei + curr_par%ech*potent_r2(r2, st%ndat)
              endif
           end if
        end if
     end do
     e = e + par%ech*ei
     pvir = pvir + par%ech*pviri
  end do
  NRG_BOX%noverlap=.true.
  NRG_BOX%ep = st%ndat%elpot_pref * e
  if(pressons) press%virp=st%ndat%elpot_pref * pvir/(3.D0*b%volume)*UnitPsT*st%inp%temp
  if(pressons) press%hsp=kT*(2*phs1-phs2)/(3.D0*b%volume)*UnitPsT*st%inp%temp
  
end function NRG_BOX

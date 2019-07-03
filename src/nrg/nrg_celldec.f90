function NRG_1P(st,b,par)
  !$ use omp_lib
  implicit none
  type(NrgRec) :: NRG_1P
  type(PcmcState),intent(in) :: st
  type(Box),intent(in),target :: b
  type(Particule),intent(in),target :: par

  integer :: cell_par,cell,i,nvois_l,counter=0,n_r
  real(8) :: e,r2,rc2
  real(8) :: prepos(3), relpos(3)
  type(Particule), pointer :: curr_par
  !$ integer :: ichunk,nthreads
  !$ nthreads = omp_get_max_threads()
  counter=counter+1
  n_r = st%dist%nfam
  rc2 = st%rcrc
  e=0._8

  cell_par = get_cell(st%clistdecomp(b%num),par%pos)
  nvois_l = st%clistdecomp(b%num)%cells(cell_par)%nvois
  !$ ichunk= nvois_l/nthreads
  !$omp parallel do default(shared),private(i,r2,curr_par,cell,prepos,relpos)&
  !$omp reduction(+:e),schedule(static,ichunk)
  do i=1,nvois_l
     cell = st%clistdecomp(b%num)%cells(cell_par)%voisin(i)
     curr_par => st%clistdecomp(b%num)%cells(cell)%hoc
     if(associated(curr_par,par)) curr_par => curr_par%suivant
     do while(associated(curr_par))
        prepos=par%pos-curr_par%pos
        r2=DIST2_MI(prepos,relpos, b%met)
        if(r2<(par%rayon+curr_par%rayon)**2) then
           NRG_1P%noverlap=.false.
#ifndef _OPENMP
           return
#endif
        elseif(r2<rc2) then
           e = e + curr_par%ech*potent_r2(r2, st%ndat)
        end if
        curr_par => curr_par%suivant
        if(associated(curr_par,par)) curr_par => curr_par%suivant
     end do
  end do
  NRG_1P%noverlap = .true.
  NRG_1P%ep = st%ndat%elpot_pref * par%ech * e

end function NRG_1P

function NRG_BOX(st,b,press)
  implicit none
  type(NrgRec) :: NRG_BOX
  type(PcmcState), intent(in) :: st
  type(Box), intent(in), target :: b
  type(Pressure),intent(out),optional :: press

  logical :: pressons
  integer :: cell_par,cell,j,i,ncell,nvois_l,n_r
  real(8) :: dump,r,ei,r2,pviri
  real(8) :: relpos(3), prepos(3),rc2, dr_press
  real(8) :: e, phs1, phs2, pvir
  type(Particule), pointer :: curr_par,par
  pressons = present(press)
  dr_press = st%inp%dr_press
  n_r = st%dist%nfam
  rc2 = st%rcrc
  e=0.D0
  phs1=0.D0
  pvir=0.D0
  phs2=0.D0
  ncell = st%clistdecomp(b%num)%ncell
  !$omp parallel do default(shared)&
  !$omp private(i,j,dump,nvois_l,r,r2,par,curr_par,cell,cell_par,prepos,relpos,ei,pviri)&
  !$omp reduction(+:e,phs1,phs2,pvir)
  do i=1,b%n
     ei = 0._8
     pviri = 0._8
     par => b%parts(i)
     cell_par = par%cellule
     nvois_l = st%clistdecomp(b%num)%cells(cell_par)%nvois
     do j=1,nvois_l
        cell = st%clistdecomp(b%num)%cells(cell_par)%voisin(j)
        curr_par => st%clistdecomp(b%num)%cells(cell)%hoc
        if(associated(curr_par,par)) curr_par => curr_par%suivant
        do while(associated(curr_par))
           prepos=par%pos-curr_par%pos
           r2=DIST2_MI(prepos, relpos,b%met)
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
           curr_par => curr_par%suivant
           if(associated(curr_par,par)) curr_par => curr_par%suivant           
        end do
     end do
     e = e + par%ech*ei
     pvir = pvir + par%ech*pviri
  end do

  NRG_BOX%noverlap = .true.
  NRG_BOX%ep = 0.5D0 * st%ndat%elpot_pref * e
  if(pressons) press%virp=0.5d0*st%ndat%elpot_pref*pvir/(3.D0*b%volume)*UnitPsT*st%inp%temp
  if(pressons) press%hsp=0.5d0*kT*(2*phs1-phs2)/(3.D0*b%volume)*UnitPsT*st%inp%temp

end function NRG_BOX

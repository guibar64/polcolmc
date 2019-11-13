module conftemp
!! Calculation related to the configurational temperature
!! T_c = k_B F/ k_B V''
use state, only: PcmcState, max_boxes
implicit none

integer, private :: step_ct_n = 0, out_ct = 38, step_ct=0

real(8), private :: force2(max_boxes), vseconde(max_boxes)

private :: potent_seconde, conf_temp1p

contains


real(8) function potent_seconde(R,i,j,st,nr)
  !! Computes the second derivative of the pair potential.
  use nrg
  implicit none
  integer,intent(in) :: i !! family
  integer,intent(in) :: j !! family
  real(8),intent(in) :: R !! separation
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  real(8), parameter :: dr = 0.01, i2dr = 1./(2*dr)
  real(8) :: kpr, kpri, kp
 
  if(st%mode_pot==1 .or. st%mode_pot==2) then ! Yukawa(2)
    kp = st%ndat%inv_debye_length
    kpr = kp*R
    kpri = 1.d0 / kpr
    potent_seconde=st%ndat%Eff_charge(i)*st%ndat%Eff_charge(j)*exp(-kpr)&
      *kp**3*kpri*(1 + 2*kpri + 2*kpri**2)
  else
    ! Could be improved.for other built_in cases
    potent_seconde=(nr%force(R+dr,i,j,st)-nr%force(R-dr,i,j,st))*i2dr 
  end if

end function


subroutine conf_temp1p(st,nr,force2,vseconde,par,n,parts,metricG,rc2)
  !! Computes, for particle `par`, the squared force `force2` and the second derivative of
  !! the potential (`vseconde`).
  use nrg
  use celldec, only: Particule
  use state, only: PcmcState
  use misc, only: dist2_mi
  implicit none
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  integer,intent(in) :: n !! number of particles
  real(8),intent(in) :: rc2 !! square of the cutoff radius
  real(8),intent(in) :: MetricG(3,3) !! box metric
  type(Particule),intent(in) :: par !! particle
  type(Particule),intent(in), dimension(n) :: parts !! all particles
  real(8),intent(out) :: force2,vseconde

  integer :: i
  real(8) :: r2,r,force2_l(3),vseconde_l,relpos(3), prepos(3)

  force2_l=0.D0
  vseconde_l=0.D0

  
  do i=1,n
     if(parts(i)%numero /= par%numero) then
     
        prepos=par%pos-parts(i)%pos
        r2=dist2_mi(prepos, relpos,MetricG,st%mode_box)
       if(r2<rc2) then
           r=sqrt(r2)
           force2_l(:)=force2_l(:)+nr%force(r,par%Famille,parts(i)%Famille, st)*relpos/r
           vseconde_l=vseconde_l+potent_seconde(r,par%Famille,parts(i)%Famille,st,nr)
        end if
     end if
  end do
  
  force2=dot_product(force2_l,force2_l)
  vseconde=vseconde_l

end subroutine conf_temp1p


subroutine conftemp_init(st, nr, fich_out_ct, first_step)
  !! initialise configurational temperature analysis
  use nrg
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  character(*), intent(in) :: fich_out_ct !! 
  integer, optional :: first_step  !! first step (default is 0)
  if(present(first_step)) then
    step_ct_n = first_step
  else
    step_ct_n = 0
  end if
  step_ct = 0
  force2 = 0.0
  vseconde = 0.0
  open(newunit=out_ct,file=fich_out_ct)

  call conftemp_new_step(st, nr, step_ct_n)
end subroutine conftemp_init

subroutine conftemp_new_step(st, nr, step)
  !! Do a new step. Appends the current quantities to the files
  use nrg
  type(PcmcState) :: st
  type(NrgRoutines) :: nr
  integer, intent(in) :: step !! step number (for printing)
  integer :: i,k
  real(8) :: force2_l(max_boxes), vseconde_l(max_boxes), fp, v2p
  step_ct = step_ct + 1
  Force2_l = 0
  Vseconde_l = 0
  do k=1, st%ntotbox
    do i=1,st%boxes(1)%n 
      call conf_temp1p(st,nr,fp,v2p,st%boxes(1)%parts(i),st%boxes(1)%n,&
           st%boxes(1)%Parts,st%boxes(1)%met,st%rcrc)
      Force2_l(k) = Force2_l(k) + fp
      Vseconde_l(k) = Vseconde_l(k) + v2p
    end do
    Force2(k) = Force2(k) + Force2_l(k)
    vseconde(k) = vseconde(k) + vseconde_l(k)
  end do
  ! Second value is in fact T_c/T
  write(out_ct,'(I10)', advance='no') step
  do k=1,st%ntotbox
    write(out_ct,'(3G16.7)', advance='no') force2_l(k)/vseconde_l(k),force2_l(k),vseconde_l(k)
  end do
  write(out_ct,*)
  
end subroutine conftemp_new_step

subroutine conftemp_finish(nb, f2, vpp)
  !! Returns average squared force and second derivative potential.
  integer, intent(in) :: nb       !! number of boxes
  real(8), intent(out) :: f2(nb)  !! squared force
  real(8), intent(out) :: vpp(nb) !! second derivative potential.
  if(step_ct > 0) then
    f2(1:nb) = force2(1:nb) / step_ct
    vpp(1:nb) = vseconde(1:nb) / step_ct
  end if
end subroutine conftemp_finish

end module conftemp

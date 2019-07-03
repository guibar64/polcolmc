module fluct_dens
!! Density and particle number fluctuations inside small volumes.
use celldec, only: Particule
use state, only: PcmcState
use agr
implicit none

logical,parameter :: use_partial_volume=.true.
real(8),parameter :: PI=3.1415926535897931D0,eps_ppho=1.D-3
integer :: N_pho,N_tests,N_min,N_N,deltaN,step_calc,N_fam,N_fam2,nf
real(8) :: Rayon,Volume,phom,pho2m,P_pho,rad_pho,phi_pho,pho_n,fluct_dens_moyf,fluct_dens_covf
allocatable :: Rayon(:),Volume(:),phom(:),P_pho(:,:),pho2m(:),rad_pho(:,:),phi_pho(:,:),N_min(:),N_N(:),&
     pho_n(:,:),deltaN(:),fluct_dens_covf(:,:),fluct_dens_moyf(:,:),N_fam(:,:),N_fam2(:,:)

private :: randomisation,dist2_rel,PI,eps_ppho,Volume,N_tests,Rayon,phom,pho2m,rad_pho,phi_pho,N_min,&
     deltaN,N_N,step_calc,P_pho,N_fam,N_fam2,nf,N_pho,fluct_dens_init_para,fluct_dens_new_step

integer, private :: fluct_dens_n_rayon !! number of sampled radii
real(8),allocatable, private :: rad_fluctu(:) !! radius fluctuations
real(8),allocatable, private :: gammas(:) !! Γ(R)

contains

subroutine fluct_dens_free()
  !! Free calculation arrays.
  implicit none
  
  if(allocated(Rayon)) then
  
    deallocate(Rayon,Volume,phom,pho2m,N_min,N_N,deltaN)

    deallocate(pho_n,P_pho,rad_pho,phi_pho)
  end if

  if(allocated(N_fam)) then
    deallocate(N_fam,N_fam2)
  end if
  
  if(allocated(fluct_dens_covf)) then
    deallocate(fluct_dens_covf,fluct_dens_moyf)
  end if

end subroutine fluct_dens_free


subroutine fluct_dens_init_para(n,L,nf_,fpara)
  !! Init calculation.
  implicit none
  integer :: n !! number of particles
  integer :: nf_ !! number of families
  real(8) :: L !! box length
  character(256) :: fpara !! parameter file
  intent(in) :: n,fpara,L

  logical :: estla,use_rel_radius
  integer :: i,j,idump,ufpara,graine
  real(8) :: pho_min,pho_max,pho_g,dump
  character(10) :: option
  
  nf=nf_

  pho_g=n/(L**3)

   ! read parameters
  inquire(file=fpara, exist=estla)
  
  if(.not.estla) then
     write(*,'(A,A)') "Error: cannot find: ", fpara
     stop 1
  endif


  open(newunit=ufpara,file=fpara,form='formatted',status='old')
  read(ufpara,*)
  read(ufpara,*) graine
  read(ufpara,*)
  read(ufpara,*) N_tests
  read(ufpara,*)
  read(ufpara,*) fluct_dens_n_rayon
  allocate(Rayon(fluct_dens_n_rayon),Volume(fluct_dens_n_rayon),phom(fluct_dens_n_rayon),pho2m(fluct_dens_n_rayon),&
       &N_min(fluct_dens_n_rayon),N_N(fluct_dens_n_rayon),deltaN(fluct_dens_n_rayon))
  read(ufpara,*)
  read(ufpara,*) option
  if(option(1:3)=="RR") then
     use_rel_radius=.true.
  elseif(option(1:2)=="R") then
     use_rel_radius=.false.
  else
     use_rel_radius=.false.
  endif
  do i=1,fluct_dens_n_rayon
     read(ufpara,*) Rayon(i)
  end do
  read(ufpara,*) 
  read(ufpara,*) pho_min
  read(ufpara,*) pho_max
  read(ufpara,*) N_pho

  close(ufpara)

  call randomisation(graine)

  do j=1,fluct_dens_n_rayon
     if(use_rel_radius) Rayon(j) = L*Rayon(j)
     Volume(j)=4.D0/3.D0*PI*(Rayon(j)**3)
  end do

  do j=1,fluct_dens_n_rayon
     dump=pho_g*Volume(j)
     N_min(j)=max(floor(dump-pho_min*sqrt(dump)),1)
     idump=floor(dump+pho_min*sqrt(dump))+1-N_min(j)
     i=idump/N_pho
     if(i<1) then
        N_N(j)=idump
        deltaN(j)=1
     else
        N_N(j)=idump/i+1
        deltaN(j) = i
     end if
  end do
  N_pho=maxval(N_N)

  allocate(pho_n(N_pho,fluct_dens_n_rayon),P_pho(n_pho,fluct_dens_n_rayon),rad_pho(n_pho,fluct_dens_n_rayon),&
&phi_pho(n_pho,fluct_dens_n_rayon))

  do j=1,fluct_dens_n_rayon
     do i=1,N_N(j)
        pho_n(i,j) = ((i-1)*deltaN(j)+N_min(j))/Volume(j)
     end do
  end do
  
  allocate(N_fam(nf,fluct_dens_n_rayon),N_fam2(nf*(nf+1)/2,fluct_dens_n_rayon))

  P_pho=0.D0
  rad_pho=0.D0
  phi_pho=0.D0
  phom(:)=0.D0
  pho2m(:)=0.D0
  N_fam=0
  N_fam2=0
  
  step_calc=0

  
end subroutine fluct_dens_init_para

subroutine fluct_dens_reset()
  !! Reset cumulative quantities. Can be used with a new run.
  implicit none

  P_pho=0.D0
  rad_pho=0.D0
  phi_pho=0.D0
  phom(:)=0.D0
  pho2m(:)=0.D0
  N_fam=0
  N_fam2=0

  step_calc=0

end subroutine fluct_dens_reset

subroutine fluct_dens_new_step(x,n,L)
  !! Performs a new step.
  implicit none
  integer,intent(in) :: n !! number of particles
  real(8),intent(in) :: L !! box length
  type(Particule),intent(in) :: x(n) !! particles
  integer :: k,j,i,idump,Ndump
  real(8) :: phi
  real(8) :: RR,dump,RR_,ray,zi,h1,h2,rad
  real(8) :: pos(3), r2
  integer :: ii,jj,Ndump_fam(nf)

  do k=1,N_tests
    do i=1,3
      call random_number(dump)
      pos(i)=L*dump
    end do
    do j=1,fluct_dens_n_rayon
      ! Nb inside the sphere
      Ndump_fam=0
      Ndump=0
      phi=0.D0
      rad=0.D0
      RR=Rayon(j)*Rayon(j)
      do i=1,n
        r2=dist2_rel(x(i)%pos,pos,L)
        if(r2<=RR) then
          Ndump=Ndump+1
          Ndump_fam(x(i)%famille)=Ndump_fam(x(i)%famille) +1
          if(.not.use_partial_volume) phi = phi + (x(i)%rayon)**3
          rad=rad+x(i)%rayon
        end if
        if(use_partial_volume) then
          dump=sqrt(r2)
          ray=x(i)%rayon
          RR_=Rayon(j)
          if(dump > RR_ .and. dump<(RR_+ray)) then
            zi=(dump*dump+RR_*RR_-ray*ray)/(2*dump)
            h2=zi+ray-dump
            h1=RR_-zi
            phi=phi+(Polyn_cal(h1,RR_)+Polyn_cal(h2,ray))
          else if(dump <= RR_ .and. dump >(RR_-ray)) then
            zi=(dump*dump+RR_*RR_-ray*ray)/(2*dump)
            h2=dump+ray-zi
            h1=RR_-zi
            phi=phi+(ray**3+Polyn_cal(h1,RR_)-Polyn_cal(h2,ray))
          else if(dump <= (RR_-ray)) then
            phi=phi+ray**3
          end if
        end if
      end do
      phi= 4.D0*PI*phi/(3.D0*Volume(j))
      rad = rad/n
      
      idump=(Ndump-N_min(j))/deltaN(j)+1
      if(idump>=1 .and. idump <= N_N(j)) then
        P_pho(idump,j)=P_pho(idump,j)+1
        phi_pho(idump,j)=phi_pho(idump,j)+phi
        rad_pho(idump,j)=rad_pho(idump,j)+rad
      end if

      phom(j)=phom(j)+Ndump
      pho2m(j)=pho2m(j)+Ndump*Ndump

      do ii=1,nf
        N_fam(ii,j)=N_fam(ii,j)*Ndump_fam(ii)
        N_fam2(ii,j)=N_fam2(ii,j)+ Ndump_fam(ii)*Ndump_fam(ii)
      end do
      idump=nf
      do ii=1,nf
        do jj=i+1,nf
          idump=idump+1
          N_fam2(idump,j)=N_fam2(idump,j)+Ndump_fam(ii)*Ndump_fam(jj)
        end do
      end do
        
    end do
    step_calc=step_calc+1

  end do
contains
  real(8) function Polyn_cal(H__,R__)
    implicit none
    real(8),intent(in) :: H__,R__
    ! Volume of a spherical cap divided by 4 PI/3
    Polyn_cal=(H__*H__*(3*R__-H__))/4
  end function Polyn_cal

end subroutine fluct_dens_new_step

subroutine fluct_dens_end_step()
  !! Finalize calculations.
  implicit none
  integer :: i,j,ind
  real(8) :: d_pho

  P_pho=P_pho/(step_calc)

  phom=phom/(step_calc)
  pho2m=pho2m/(step_calc)
  pho2m=(pho2m-phom*phom)
 
  phi_pho(:,:)=phi_pho(:,:)/step_calc
  rad_pho(:,:)=rad_pho(:,:)/step_calc
  
  where(P_pho>EPS_ppho)
    rad_pho= rad_pho/P_pho
  end where

  where(P_pho>EPS_ppho)
    phi_pho= phi_pho/P_pho
  end where

  ! Normalization
  do i=1,fluct_dens_n_Rayon
    d_pho = pho_n(2,i)-pho_n(1,i)
    P_pho(:,i)=P_pho(:,i)/d_pho
  end do

  allocate(fluct_dens_covf((nf*(nf+1))/2,fluct_dens_n_rayon),fluct_dens_moyf(nf,fluct_dens_n_rayon))

  fluct_dens_moyf(1:nf,:)=(1.D0*N_fam(1:nf,:))/step_calc
  fluct_dens_covf(1:nf,:)=(1.D0*N_fam2(1:nf,:))/step_calc-fluct_dens_moyf(1:nf,:)**2

  ind=0
  do i=1,nf
    do j=i+1,nf
      ind=ind+1
      fluct_dens_covf(ind,:)=(1.D0*N_fam2(ind,:)/step_calc)-(fluct_dens_moyf(i,:)*fluct_dens_covf(j,:))
    end do
  end do


  
end subroutine fluct_dens_end_step

subroutine fluct_dens_get_gammas(gams)
  !! Get 
  implicit none
  real(8), intent(out) :: gams(*)
  gams(1:fluct_dens_n_rayon)=pho2m(1:fluct_dens_n_rayon)/phom(1:fluct_dens_n_rayon)
end subroutine fluct_dens_get_gammas

subroutine fluct_dens_get_radii(rads)
  implicit none
  real(8), intent(out) :: rads(*)
  rads(1:fluct_dens_n_rayon)=rayon(1:fluct_dens_n_rayon)
end subroutine fluct_dens_get_radii

subroutine randomisation(seed_aucasou)
  !! Initialize fortran's random generator.
  implicit none
  integer,intent(in) :: seed_aucasou
  
  integer,parameter :: MX=2047483647
  integer :: size_seed,seed,i,idump
  real(8) :: dump
  allocatable :: seed(:)

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  call system_clock(count=idump,count_rate=idump)

  if(idump<=0) then
     idump=seed_aucasou
   end if
   seed=idump
   call random_seed(put=seed)
   do i=2,size_seed
     call random_number(dump)
     seed(i)=floor(MX*dump)
   end do
   call random_seed(put=seed)
  

   deallocate(seed)
end subroutine randomisation

real(8) function dist2_rel(coord0,coord,L)
  !! square of the separation betwen two images
  ! TODO: use the real dist_mi ? 
  implicit none
  real(8),intent(in) :: coord0(3),coord(3),L
  real(8) :: dump,Lsur2
  integer :: i
  Lsur2=0.5D0*L
  
  dist2_rel=0
  do i=1,3
    dump=coord(i)-coord0(i)
    if(dump>Lsur2) then
      dump = dump-L
    elseif(dump<-lsur2) then
      dump = dump+L
    end if
    dist2_rel = dist2_rel+dump*dump
  end do

end function dist2_rel

subroutine fluct_dens_init(st, file_name)
  !! Initialize calculations
  type(PcmcState) :: st
  character(*), intent(in) :: file_name !! configuration file
  logical :: inestla
  inquire(file=st%inp%fich_fluctu,exist=inestla)
  if(.not.inestla) then
    write(*,'(A,A,A)') "Error", trim(file_name)," not found."
    stop 1
  endif

  call fluct_dens_init_para(st%boxes(1)%n,real(st%boxes(1)%tens(1,1),8),st%dist%nfam,st%inp%fich_fluctu)
end subroutine fluct_dens_init

subroutine fluct_dens_update(st)
  !! Does a new calculation step from current state.
  type(PcmcState) :: st
  integer :: k
  do k=1,1
    call fluct_dens_new_step(st%boxes(1)%Parts,st%boxes(1)%n,st%boxes(1)%tens(1,1))
  end do
end subroutine fluct_dens_update

subroutine fluct_dens_finish()
  !! Finalize calculations
  call fluct_dens_end_step()
  
  allocate(gammas(fluct_dens_n_rayon),rad_fluctu(fluct_dens_n_rayon))
  call fluct_dens_get_gammas(gammas)
  call fluct_dens_get_radii(rad_fluctu)
end subroutine fluct_dens_finish

subroutine fluct_dens_output(fdist,frad,fphi,fgg)
  !! Prints results to files
  character(*),intent(in) :: fdist !! Density distributions
  character(*),intent(in) :: frad  !! Mean radii
  character(*),intent(in) :: fphi  !! Volume fractions
  character(*),intent(in) :: fgg   !! Fluctuations

  integer :: i,j,k,idump,ufout
  
  open(newunit=ufout,file=fdist,form="formatted")
  call agr_init(ufout)
  write(ufout,'(A)') '@xaxis  label "\\xr"'
  write(ufout,'(A)') '@yaxis  label "Dist"'
  do i=1,fluct_dens_n_rayon
    write(ufout,'(A,I0,A,G15.7,A)') "@s",i-1,' legend "',Rayon(i),' nm"'
  end do
  do i=1,fluct_dens_n_rayon
    call agr_add_xy(ufout,i-1,pho_n(:,i),P_pho(:,i),N_N(i))
  end do
  close(ufout)

  open(newunit=ufout,file=frad,form="formatted")
  call agr_init(ufout)
  write(ufout,'(A)') '@xaxis  label "\\xr"'
  write(ufout,'(A)') '@yaxis  label "Mean radius"'
  do i=1,fluct_dens_n_rayon
    write(ufout,'(A,I0,A,G15.7,A)') "@s",i-1,' legend "',Rayon(i),' nm"'
  end do
  do i=1,fluct_dens_n_rayon
    call agr_add_xy(ufout,i-1,pho_n(:,i),rad_pho(:,i),N_N(i))
  end do
  close(ufout)

  open(newunit=ufout,file=fphi,form="formatted")
  call agr_init(ufout)
  write(ufout,'(A)') '@xaxis  label "\\xr"'
  write(ufout,'(A)') '@yaxis  label "mean \\xf\\f{}"'
  do i=1,fluct_dens_n_rayon
    write(ufout,'(A,I0,A,G15.7,A)') "@s",i-1,' legend "',Rayon(i),' nm"'
  end do
  do i=1,fluct_dens_n_rayon
    call agr_add_xy(ufout,i-1,pho_n(:,i),phi_pho(:,i),N_N(i))
  end do
  close(ufout)
  

  open(newunit=ufout,file=fgg,form="formatted")

  write(ufout,*) "Fluctuations of families: cov(Ni,Nj)"
  write(ufout,*) "Nb radii:",fluct_dens_n_rayon
  write(ufout,*) "Nb families:",nf
  do k=1,fluct_dens_n_rayon
    write(ufout,*) Rayon(k)
    idump=0
    do i=1,nf
      idump=idump+1
      write(ufout,'(I16,1x,G15.7)') i,fluct_dens_covf(idump,k)
    end do
    do i=1,nf
      do j=i+1,nf
        idump=idump+1
        write(ufout,'(I16,1x,I16,1x,G15.7)') i,j,fluct_dens_covf(idump,k)     
      end do
    end do
  end do

  close(ufout)

end subroutine fluct_dens_output

subroutine fluct_dens_log_results(st, log_unit)
  !! Prints main results.
  type(PcmcState) :: st
  integer, intent(in) :: log_unit !! file descriptor
  integer :: i,k
  
  write(log_unit,*) "Density Fluctuations:"
  do k=1, st%ntotbox
    write(log_unit,*) "Box", k
    write(log_unit,*) "phi   Rayon   Gamma^-1"
    do i=1,fluct_dens_n_rayon
      write(log_unit,'(3G16.8)') st%res(k)%hetv,rad_fluctu(i),gammas(i)
    end do
  end do
  
end subroutine fluct_dens_log_results


end module fluct_dens

subroutine CORRELS_NEW_STEP_DEF(corrs,metricG,n,centres)
  implicit none
  type(CorrelsRec),intent(inout) :: corrs
  integer,intent(in) :: n
  real(8),intent(in),dimension(3,3) :: metricG
  type(Particule),intent(in) :: centres(*)

  integer :: i,j,k,indexg,ng,nr
  real(8) :: r2,prepos(3),relpos(3),rmin,dr, dri

  nr=corrs%nr
  ng=corrs%ng
  rmin=corrs%rmin
  dr = corrs%dr
  dri = 1.d0/dr

  !$omp parallel do default(shared) private(i,j,prepos, relpos, indexg)
  do i=1,n
     do j=i+1,n
        prepos=centres(i)%pos-centres(j)%pos
        r2=sqrt(DIST2_MI(prepos,relpos,metricG))
        
        k=Indic_(centres(i)%famille,centres(j)%famille,nr)
        indexg=1+floor((r2-rmin)*dri)
        if(indexg>0.and.indexg<=ng) then
          !$omp atomic update
          corrs%nder(indexg,k)=corrs%nder(indexg,k)+2
        endif
     end do
  end do
  !$omp end parallel do
  corrs%step_calc=corrs%step_calc+1

end subroutine


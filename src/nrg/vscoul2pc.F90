module vscoul_mod
  !! Fast (and dirty) implementation for r² → exp(r²)/r²
  !!
  !! This relie on IEEE float trickery, the steps are approximately:
  !!```
  !!x → ≈ log₂ x + cte → index → c(index)
  !!```
  !! where c is a precomputed table. 
  use gaussmod
  implicit none
  integer nPol, nShift, MxTab
  parameter (nPol=4, nShift=16)
  parameter (MxTab=2000)

  real(8) :: c(0:MxTab)

#if POLCOLMC_BIG_ENDIAN == 1
  integer, parameter :: imsw=1,ilsw=2
#else
  integer, parameter :: imsw=2,ilsw=1
#endif
  integer :: indx,ind0,ind1,ind2
  integer :: IndMin,IndMax
  data indx/z'3f170000'/, ind0/z'3f180000'/ ,ind1/z'408F4000'/ , ind2/z'00010000'/       

  private :: fnc, MkErr, MkFit

contains
  function fnc(x)
    !! The function to be tabulated.
    real(8) :: fnc
    real(8), intent(in) :: x
    fnc=exp(-sqrt(x))/sqrt(x)
  end function fnc

  subroutine GenTab(xmin,xmax, lout)
    !! Generates the table. Sets `xmin` and `xmax`.
    real(8), intent(out) :: xmin !! Minimum accessible by the table 
    real(8), intent(out) :: xmax !! Maximum accessible by the table 
    logical, intent(in) :: lout  !! If set to `.true.`, some information about the table is printed on standard output.
    integer :: ix0, ix1, ixx
    real(8) :: x0, x1, ErrMax, err
    integer :: i, index
    dimension :: ix0(1:2),ix1(1:2),ixx(1:2)

    ix0(imsw)=ind0
    ix0(ilsw)=0
    ix1(imsw)=ind1
    ix1(ilsw)=0
    x0 = transfer(ix0, 1.d0)
    x1 = transfer(ix1, 1.d0)
    xmin = x0
    xmax = x1
    index=0
    do i=0,nPol
      c(i)=0.0d0
    end do
    c(1) =-fnc(x0)*(1.0d0/sqrt(x0)+1.0d0/x0)
    c(0) = fnc(x0)-c(1)*x0
    index=index+nPol+1
    ErrMax=0.0d0
    if(lout) then
      write(*,'(a,/)')'********** VSCIUL2PC.F90: GenTab routine ***********'
      Write(*,'(a,f12.8)') ' Smallest argument:         ',x0
      Write(*,'(a,f12.8)') ' Its value:                 ',fnc(x0)
      Write(*,'(a,g15.8)') ' Largest argument:          ',x1
      Write(*,'(a,f12.8)') ' Its value:                 ',fnc(x1)
      Write(*,'(a,i12)') ' Polynomial degree:         ',nPol
      write(*,*)
    endif
    do i=ind0,ind1,ind2
      ix0(imsw)=i
      ix1(imsw)=i+ind2
      x0 = transfer(ix0, 1.d0)
      x1 = transfer(ix1, 1.d0)
      call MkFit(x0,x1,c(index),nPol)
      call MkErr(x0,x1,c(index),nPol,err)
      ErrMax=max(ErrMax,err)
      index=index+nPol+1
    end do
    index=index-nPol-1
    do i=0,nPol
      c(i+index)=0.0d0
    end do
    index=index+nPol+1
    if(lout) then
      write(*,'(a,i12)') ' Size of table:             ',index 
      write(*,'(a,e10.3,/)') ' Max error of algorithm:      ',ErrMax
    end if
    IndMax=iShft(ind1-indx,-nShift)
    IndMin=0
  end subroutine GenTab

  subroutine MkErr(x0,x1,Coef,nPol,ErrMax)
    real(8) :: x0, x1, Coef, ErrMax
    integer :: nPol
    dimension :: Coef(0:nPol)
    real(8) :: Err, Exact, Fit, t, x
    integer :: i, ifin, k
    ErrMax=0.0d0
    ifin = int(1/(0.01d0)+1)
    do i=1,ifin,1
      x=x0+(i-1)*0.01d0*(x1-x0)
      Exact=fnc(x)
      t=0.0d0
      do k=nPol,0,-1
        t=t*x+Coef(k)
      end do
      Fit=t
      Err=Exact-Fit
      ErrMax=max(ErrMax,abs(Err))
    end do
  end subroutine MkErr

  subroutine MkFit(x0,x1,Coef,nPol)
    real(8) :: x0,x1,Coef
    integer :: nPol
    integer :: MxPol
    real(8) :: Hess, Absc, Fit, Exact, pi, t
    integer :: i,j,k
    parameter (MxPol=5)
    dimension  :: Coef(0:nPol)
    dimension  :: Hess(0:MxPol,0:MxPol)
    dimension :: Absc(0:MxPol)
    dimension :: Fit(0:MxPol),Exact(0:MxPol)
    pi=4.0d0*atan(1.0d0)
    do i=0,nPol
      Absc(i)=0.5d0*(x0+x1)-0.5d0*(x1-x0)*cos(pi*(2*i+1)/(2*nPol+2))
    end do
    do i=0,nPol
      t=1.0d0
      do j=0,nPol
        Hess(i,j)=t
        t=t*Absc(i)
      end do
      Coef(i)=fnc(Absc(i))
    end do
    call Gauss(nPol+1,MxPol+1,Hess,Coef,Coef)
    do i=0,nPol
      Exact(i)=fnc(Absc(i))
      t=0.0d0
      do k=nPol,0,-1
        t=t*Absc(i)+Coef(k)
      end do
      Fit(i)=t
    end do
  end subroutine MkFit


end module vscoul_mod

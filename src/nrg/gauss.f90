module gaussmod
!! Helper for [[vscoul_mod]]
contains
  
subroutine gauss(n,ldim,a,x,c)
   implicit none
   integer n,ldim
   real(8) a,x,c
   integer i,j,k
   real(8) swap,fact
   dimension a(ldim,ldim),x(ldim),c(ldim)
   do i=1,n
      x(i)=c(i)
   end do
   do i=1,n-1
      k=i
      do j=i+1,n
         if(abs(a(k,i)) .lt. abs(a(j,i))) k=j
      end do
      if(k .ne. i) then
         do j=i,n
            swap=a(i,j)
            a(i,j)=a(k,j)
            a(k,j)=swap
         enddo
         swap=x(i)
         x(i)=x(k)
         x(k)=swap
      end if
      do k=i+1,n
         fact=a(k,i)/a(i,i)
         do j=i+1,n
            a(k,j)=a(k,j)-fact*a(i,j)
         end do
         x(k)=x(k)-fact*x(i)
      end do
   end do
   x(n)=x(n)/a(n,n)
   do i=n-1,1,-1
      do k=i+1,n
         x(i)=x(i)-a(i,k)*x(k)
      end do
      x(i)=x(i)/a(i,i)
   end do
   return
end subroutine

end module gaussmod

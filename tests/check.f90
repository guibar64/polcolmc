logical function appr(x,y)
real(8) :: x,y
if(abs(x-y) < 1.d-6*abs(x)) then
    appr = .true.
else
    appr = .false.
end if
end function

subroutine check_double(x,y)
real(8) :: x,y
if(.not. appr(x,y)) then
    write(*,'(G16.7,"/=",G16.7)') x,y
    stop 1
end if
end subroutine
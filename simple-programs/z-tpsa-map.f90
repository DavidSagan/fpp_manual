module small_code 
use pointer_lattice

implicit none
 INTERFACE T
     MODULE PROCEDURE track_real
     MODULE PROCEDURE track_real_8
 END INTERFACE
contains 

subroutine track_real(xs0)
implicit none
type(probe) xs0
 xs0%x(1)=xs0%x(1)+1.5_dp*xs0%x(1)
 xs0%x(2)=xs0%x(2)-0.3_dp*xs0%x(1) + sin(xs0%x(1))**2
end subroutine track_real

subroutine track_real_8(xs0)
implicit none
type(probe_8) xs0
 xs0%x(1)=xs0%x(1)+1.5_dp*xs0%x(1)
 xs0%x(2)=xs0%x(2)-0.3_dp*xs0%x(1) + sin(xs0%x(1))**2
end subroutine track_real_8

end module small_code 

program example1
use madx_ptc_module
use pointer_lattice
use small_code
implicit none
type(probe) xs0
type(probe_8) xs
real(dp) prec,x(6)
integer i
longprint=.false.
 
call ptc_ini_no_append
call append_empty_layout(m_u)
call set_pointers; use_quaternion=.true.;
x=0.001d0
xs0=x

call T(xs0)
call print(xs0)

call alloc(xs)
x=0.001d0
xs0=x
xs=xs0
call T(xs)
call print(xs)

call init(only_2d0,1,0) ! TPSA set no=1 and nv=2

x=0.0d0
xs=x

write(6,*) " Must be real "
call print(xs%x(1))
do i=1,2
 xs%x(i)=xs0%x(i)+dz_8(i)
enddo
write(6,*) " Must be Taylor "
call print(xs%x(1))

call T(xs)
write(6,*) " Must part of a potential map for the system "
call print(xs%x(1))
call print(xs)

end program example1








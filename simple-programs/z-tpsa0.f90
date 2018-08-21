program example0
use madx_ptc_module
use pointer_lattice
implicit none
type(real_8) x_8
real(dp) x
longprint=.false.
call ptc_ini_no_append
call append_empty_layout(m_u)
call set_pointers; use_quaternion=.true.;

call init(only_2d0,3,0) ! TPSA set no=3 and nv=2

x=0.1d0

call alloc(x_8)
 
x_8=x  ! insert a real into a polymorph
 
write(6,*) " Must be real "
call print(x_8)
 
x_8=x_8+dz_8(1)
 
write(6,*) " Must be Taylor : 0.1+dx "
call print(x_8)

x_8=x_8**4
 
write(6,*) " Must (0.1+dx)**4 to third order "
 
call print(x_8)
end program example0








program my_small_code_real_8
use my_code
use c_tpsa
implicit none
real(dp) :: z0(2) = (/0,0/)  ! special orbit
type(probe) :: probe0 
type(probe_8) :: polymorphic_probe  
type(c_damap) Identity, one_period_map,diag,two_turn_map,one_period_map_squared
type(c_normal_form) normal
integer :: nd = 1 , no = 2 , np = 0 , imap , ipol,again
logical success,ok
longprint=.false.
! nd = number of degrees of freedom
! no =  order of Taylor series
! Number of extra variables beyond 2*nd

call read_lattice ("lattice-2.txt")

100 continue
print * , " Do you want polymorphs as parameters of the Taylor series ?"
print * , " 0=> no, 1=> yes ?"; read(5,*) ipol
if(ipol==1) call read_polymorphs 
call make_polymorphs(np) !
call find_closed_orbit(z0,ok) 
if(ok) then
 print * ,"closed orbit found"
else
 print * ,"closed orbit NOT found"
stop
endif

print * ," Closed Orbit",z0
print * , "The order of the Taylor series ?"
read(5,*) no

call c_init_all(no,nd,np) ! Initializes TPSA 
call alloc(polymorphic_probe)
call alloc(Identity, one_period_map,diag,two_turn_map, &
one_period_map_squared) ; call alloc(normal);

probe0=z0
Identity=1

polymorphic_probe = probe0 + Identity
 

call track(polymorphic_probe)

one_period_map=polymorphic_probe

call track(polymorphic_probe)

two_turn_map=polymorphic_probe

one_period_map_squared=one_period_map**(2)

Identity=one_period_map_squared**(-1)*two_turn_map
call clean(Identity,Identity,prec=1.d-9)
 call print(Identity)

print * , "Do you want to redo this calculation?"
print * , " 0=> no, 1=> yes ?"; read(5,*) again
if(again==1) goto 100

end program my_small_code_real_8



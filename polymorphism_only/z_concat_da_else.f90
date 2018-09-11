program my_small_code_real_8
use my_code
use c_tpsa
implicit none
real(dp) :: z0(2) = (/0,0/),z1(2) = (/0,0/)  ! special orbit
type(probe) :: probe0 
type(probe_8) :: polymorphic_probe  
type(c_damap) Identity,one_period_map,newton_map,to_closed_orbit
type(c_damap) two_period_map
type(c_vector_field) f
integer :: nd = 1 , no = 2 , np = 0 , imap , ipol,again
logical :: ok=.true.
longprint=.false.
! nd = number of degrees of freedom
! no =  order of Taylor series
! Number of extra variables beyond 2*nd

call read_lattice ("lattice-2.txt")

100 continue
!print * , " Do you want polymorphs as parameters of the Taylor series ?"
!print * , " 0=> no, 1=> yes ?"; read(5,*) ipol
!if(ipol==1) call read_polymorphs 
call make_polymorphs(np) !
call find_closed_orbit(z0,ok) 



!print * ," Closed Orbit from Tracking",z0
1 print * , "The order of the Taylor series ?"
read(5,*) no
print * , "Tpsa =0 DA=1 "
read(5,*) ipol
call c_init_all(no,nd,np) ! Initializes TPSA 
call alloc(polymorphic_probe)
call alloc(Identity,one_period_map,newton_map,newton_map,to_closed_orbit)
call alloc(two_period_map)
call alloc(f)
probe0=z0
Identity=1

polymorphic_probe = probe0 + Identity

 
call print(one_period_map)

call track(polymorphic_probe)



one_period_map=polymorphic_probe
write(6,*) "more"
read(5,*) again
if(again==1) goto 1
end program my_small_code_real_8



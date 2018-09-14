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

if(ok) then
 print * ,"closed orbit found"
else
 print * ,"closed orbit NOT found"
stop
endif

print * ," Closed Orbit from Tracking",z0
print * , "The order of the Taylor series ?"
read(5,*) no
print * , "Tpsa =0 DA=1 "
read(5,*) ipol
if(ipol==0) then
 z1=z0
 z0=0.0001d0
endif
call c_init_all(no,nd,np) ! Initializes TPSA 
call alloc(polymorphic_probe)
call alloc(Identity,one_period_map,newton_map,newton_map,to_closed_orbit)
call alloc(two_period_map)
probe0=z0
Identity=1

polymorphic_probe = probe0 + Identity
 

call track(polymorphic_probe)

one_period_map=polymorphic_probe

if(ipol==1) then
call PRINT(one_period_map)
elseif(ipol==2) then
 two_period_map=one_period_map*one_period_map
 print *, " one_period_map*one_period_map "
 call PRINT(two_period_map)
 call track(polymorphic_probe)
 two_period_map=polymorphic_probe
 print *, " Two periods map "
 call PRINT(two_period_map)
stop
else
newton_map=one_period_map; newton_map=maketpsa(newton_map)
newton_map%v(1)=newton_map%v(1)-(1.0_dp.cmono.1)
newton_map%v(2)=newton_map%v(2)-(1.0_dp.cmono.2)

newton_map=newton_map.oo.(-1)
print * ," Closed Orbit from Tracking ",z1
z0(1)=newton_map%v(1)
z0(2)=newton_map%v(2)
 print * ," Closed Orbit using TPSA map",z0 

to_closed_orbit%v(1)=identity%v(1)+z0(1)
to_closed_orbit%v(2)=identity%v(2)+z0(2)
 
 print * ," Linear using TPSA map around z=(0,0)"
newton_map=one_period_map.cut.2
call print(newton_map)
one_period_map= one_period_map.o.to_closed_orbit
one_period_map=(to_closed_orbit.oo.(-1)).o.one_period_map
newton_map=one_period_map.cut.2

 print * ," Linear map around computed by TPSA inversion " 

call print(newton_map)

 endif

!call clean(one_period_map,one_period_map,prec=1.d-9)
! call print(one_period_map)
call alloc(f)
!f=log(one_period_map)

!one_period_map=maketpsa(one_period_map)

f=c_logf_spin( one_period_map)  !,h,epso,n,tpsa )     !log(one_period_map)

one_period_map=exp(f,one_period_map)
call PRINT(f)
call PRINT(one_period_map)

print * , "Do you want to redo this calculation?"
print * , " 0=> no, 1=> yes ?"; read(5,*) again
if(again==1) goto 100

end program my_small_code_real_8



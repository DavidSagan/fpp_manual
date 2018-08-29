program my_small_code_real_8
use my_code
use c_tpsa
implicit none
real(dp) :: z0(2) = (/0,0/)  ! special orbit
type(probe) :: probe0 
type(probe_8) :: polymorphic_probe  
type(c_damap) Identity, one_period_map,diag
type(c_normal_form) normal
integer :: nd = 1 , no = 2 , np = 0 , imap , ipol,again
logical success,ok
! nd = number of degrees of freedom
! no =  order of Taylor series
! Number of extra variables beyond 2*nd

call read_lattice ("lattice-1.txt")

100 continue
print * , " Do you want polymorphs as parameters of the Taylor series ?"
print * , " 0=> no, 1=> yes ?"; read(5,*) ipol
if(ipol==1) call read_polymorphs 
call make_polymorphs(np) !
call find_closed_orbit(z0,ok) 
print * , "The order of the Taylor series ?"
read(5,*) no

call c_init_all(no,nd,np) ! Initializes TPSA 
call alloc(polymorphic_probe)
call alloc(Identity, one_period_map,diag) ; call alloc(normal);

probe0=z0
Identity=1

polymorphic_probe = probe0 + Identity
 

call track(polymorphic_probe)

one_period_map=polymorphic_probe
call c_normal(one_period_map,normal)
call print(one_period_map)
call print(normal%atot)
diag=c_phasor(-1)*normal%atot**(-1)*one_period_map*normal%atot*c_phasor()
call clean(diag, diag ,prec=1.e-10_dp)
call print(diag)
write(6,*)  normal%tune(1)*twopi
print * , "Do you want to redo this calculation?"
print * , " 0=> no, 1=> yes ?"; read(5,*) again
if(again==1) goto 100

end program my_small_code_real_8



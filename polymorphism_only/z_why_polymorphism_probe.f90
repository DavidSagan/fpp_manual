program my_small_code_real_8
use my_code
implicit none
type(probe) :: real_probe  
type(probe_8) :: polymorphic_probe  
real(dp) :: z0(2) = (/0,0/)  ! special orbit
integer :: nd = 1 , no = 2 , np = 0 , imap , ipol,again
logical success
! nd = number of degrees of freedom
! no =  order of Taylor series
! Number of extra variables beyond 2*nd

call read_lattice ("lattice-1.txt")

100 continue
print * , " Do you want polymorphs as parameters of the Taylor series ?"
print * , " 0=> no, 1=> yes ?"; read(5,*) ipol
if(ipol==1) call read_polymorphs 
call make_polymorphs(np) !

print * , "The order of the Taylor series ?"
read(5,*) no

call init(no,nd,np) ! Initializes TPSA 
call alloc(polymorphic_probe)

print * , " Do you want bona fide Taylor approximation of the map ?"
print * , " 0=> no, 1=> yes ?"; read(5,*) imap

polymorphic_probe%x(1)=z0(1) + imap*dz_8(1)
polymorphic_probe%x(2)=z0(2) + imap*dz_8(2)

call track(polymorphic_probe)

call print(polymorphic_probe)

print * , "Do you want to redo this calculation?"
print * , " 0=> no, 1=> yes ?"; read(5,*) again
if(again==1) goto 100

end program my_small_code_real_8
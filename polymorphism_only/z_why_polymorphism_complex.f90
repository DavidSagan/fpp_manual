program my_small_code_complex_8
use polymorphic_complextaylor
implicit none
type(complex_8) :: z(2)  
real(dp) :: z0(2) = (/0,0/)  ! special orbit
type(complex_8)  :: L  , K_s 
type(real_8)  ::   K_q  

integer :: nd = 1 , no = 2 , np = 0 ,ip

! nd = number of degrees of freedom
! no =  order of Taylor series
! Number of extra variables beyond 2*nd

call alloc(z)
call alloc(  K_q  )
call alloc( L  , K_s )
np=0
print * , "Give  L and parameter ordinality (0 if not a parameter)"
read(5,*) L%r , ip
np=np+ip
call make_it_knob(L,ip,0);  np=np+ip;
print * , "Give  K_q and parameter ordinality (0 if not a parameter)"
read(5,*) K_q%r , ip
call make_it_knob(K_q,ip);  np=np+ip;
print * , "Give  K_s and parameter ordinality (0 if not a parameter)"
read(5,*) K_s%r , ip
call make_it_knob(K_s,ip,0); np=np+ip;
print * , "The order of the Taylor series ?"
read(5,*) no



call init(no,nd,np) ! Initializes TPSA 
call print(k_q)
z(1)=k_q
call print(z(1))
stop
z(1)=z0(1) + dz_8(1)
z(2)=z0(2) + dz_8(2)

call track(z)

call print(z)

contains

subroutine track(z)
implicit none
type(complex_8) :: z(2) 
 z(1)=z(1)+L*z(2) 
 z(2)=z(2)-K_q*z(1)-K_s*z(1)**2 
end subroutine track

end program my_small_code_complex_8
      module my_code
      use tree_element_module
      implicit none
      private  trackr,trackp
      type(real_8)  :: L ,B, K_q , K_s 
      real(dp) :: L0 , B0, K_q0 , K_s0 
      real(dp) par(4)
      integer ip(4)
      
      interface track
       module procedure trackr
       module procedure trackp 
      end interface
      
      contains
       
subroutine read_lattice (filename)
implicit none
character(*), optional :: filename
integer mf,ier
call set_da_pointers()
  mf=5 ; ip=0;
 if(present(filename)) then 
  mf=99
  open(unit=mf, file=filename, status='old',iostat=ier)
  if(ier/=0) then
    print * , "No lattice file !"
    stop
  endif
  
 endif
 if(mf==99) then
  read(mf,*) par(1) 
  read(mf,*) par(2) 
  read(mf,*) par(3) 
  read(mf,*) par(4) 
  print * , "L"
  print * ,  par(1) 
  print * , "K_b"
  print * ,  par(2) 
  print * , "K_q"
  print * ,  par(3)  
  print * , "K_s"
  print * ,  par(4)  
  print * , "Lattice Read and polymorphs read "
 close(mf)
 else
  print * , "Give L"
  read(mf,*) par(1) 
  print * , "Give  K_b"
  read(mf,*) par(2) 
  print * , "Give  K_q"
  read(mf,*) par(3)  
  print * , "Give  K_s"
  read(mf,*) par(4)  
 endif

end subroutine read_lattice


subroutine read_polymorphs (filename)
implicit none
character(*), optional :: filename
integer mf,ier
  mf=5 ; ip=0;
 if(present(filename)) then 
  mf=99
  open(unit=mf, file=filename, status='old',iostat=ier)
  if(ier/=0) then
    print * , "No lattice file !"
    stop
  endif
 endif
  print * , "Give parameter ordinality of L  (0 if not a parameter)"
  read(mf,*)  ip(1)
  print * , "Give parameter ordinality of K_b " 
  read(mf,*)  ip(2)
  print * , "Give parameter ordinality of K_q "
  read(mf,*)  ip(3)
  print * , "Give parameter ordinality of K_s"
  read(mf,*)  ip(4)
end subroutine read_polymorphs


subroutine make_polymorphs(np)
implicit none
 integer   k,np
 
 call alloc( L , B, K_q , K_s )
 use_quaternion=.true.
 np=0
 L0=par(1)
 B0=par(2)
 K_q0=par(3)
 K_s0=par(4)
 L=L0
 call make_it_knob(L,ip(1));
 B=b0
 call make_it_knob(L,ip(2));  
 K_q=K_q0
 call make_it_knob(K_q,ip(3));  
 K_s=K_s0
 call make_it_knob(K_s,ip(4)); 

 do k=1,4
  if(ip(k)>np) np=k
 enddo

end subroutine make_polymorphs


      subroutine trackr(p) ! for probe
      implicit none
      type(probe) :: p
       p%x(1)=p%x(1)+L0*p%x(2) 
       p%x(2)=p%x(2)-B0-K_q0*p%x(1)-K_s0*p%x(1)**2 
      end subroutine trackr
      
      subroutine trackp(p) ! for probe_8
      implicit none
      type(probe_8) :: p
       p%x(1)=p%x(1)+L*p%x(2) 
       p%x(2)=p%x(2)-B-K_q*p%x(1)-K_s*p%x(1)**2 
      end subroutine trackp

subroutine find_closed_orbit(z0,ok) 
implicit none
real(dp) z0(2)
type(probe) :: pf,pm1,pm2,p
real(dp) :: eps=1e-6_dp,eps0=1e-7_dp
real(dp) mat(2,2),matinv(2,2),det,v(2),norm,normb
integer i,k
logical ok
ok=.true.
pf=0
pf%x(1:2)=z0
pm1=0
pm2=0
p=0
mat=0
k=0
normb=1.d38
 do i=1,100
   p=pf
   pm1=pf
   pm2=pf
   pm1%x(1)=pm1%x(1)+eps
   pm2%x(2)=pm2%x(2)+eps
   call track(p)
   call track(pm1)
   call track(pm2)
   mat(1,1)=1.0_dp-(pm1%x(1)-p%x(1))/eps
   mat(1,2)=-(pm2%x(1)-p%x(1))/eps
   mat(2,1)=-(pm1%x(2)-p%x(2))/eps
   mat(2,2)=1.0_dp-(pm2%x(2)-p%x(2))/eps
   det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
   matinv(1,1)=mat(2,2)/det
   matinv(2,2)=mat(1,1)/det
   matinv(1,2)=-mat(1,2)/det
   matinv(2,1)=-mat(2,1)/det
   v=p%x(1:2)-pf%x(1:2)
   v=matmul(matinv,v)
   pf%x(1:2)=pf%x(1:2)+v
   norm=abs(v(1))+abs(v(2)) 
   if((norm > eps0 .and. k==0) .or. i<3 ) then
     normb=norm
   elseif(norm<eps0.and.k==0) then
     k=1
     normb=norm
   elseif(k==1) then
     if(normb>=norm) exit
     normb=norm
   endif
 enddo

if(i>=100) then
 ok=.false.
else
 z0=pf%x(1:2)
endif
end subroutine find_closed_orbit 

end module my_code














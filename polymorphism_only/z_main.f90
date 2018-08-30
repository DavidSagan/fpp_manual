!  Console1.f90 
!
!  FUNCTIONS:
!  Console1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Console1
use polymorphic_complextaylor
    implicit none
    type(real_8) x
    type(complextaylor) z1,z2
     longprint=.false.
    call init(2,4)
   
    call alloc(x)
    call alloc(z1,z2)

    x=(1.d0+(1.d0.mono.1))**2*(1.d0.mono.2)

!call PRINT(x)

    z1=dz_t(1)+i_*dz_t(2)

    z2=z1**2
   write(6,*) " this is z1 "
   call PRINT(z1)   
   write(6,*) " this is z2 "
   call PRINT(z2)
    end program Console1


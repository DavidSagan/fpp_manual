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

    program tpsa_analysis
!    use tree_element_MODULE
    use c_TPSA
    implicit none
    type(real_8) x
    type(c_damap) c_map
    integer i,nd,no

    nd=2 ! degrees of freedom
    no=3 ! order of TPSA
    call c_init_all(no,nd,np1=0,ndpt1=0,AC_RF=0) 
 
    call alloc(x); call alloc(c_map)

    c_map=1 ! map in made to identity

    x=(1.d0+(1.d0.mono.1))**2*(1.d0.mono.2)

     call PRINT(c_map)
     c_map=0
     call PRINT(c_map)
     do i=1,4
       c_map%v(i)=dz_c(i)
     enddo
     call PRINT(c_map)
    end program tpsa_analysis


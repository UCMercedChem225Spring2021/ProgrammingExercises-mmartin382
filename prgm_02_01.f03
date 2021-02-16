!
        Program prgm_02_01
!
!	This program evaluates and prints kinetic energy matrix 
!	elements for particle-in-a-box (1D) eigenstates. 
!
!	Because we are looking at the PIB with potential 
!	energy V=0, the total energy is only comprised of T.
!	Additionally, orthonormality requires that we will
!	only obtain a nonzero value when n1=n2.
!
!	Atomic units means h=1.
!
        implicit none
        real :: m,l
        integer :: n1,n2
        integer :: ierror = 1
	real :: PIB_1D_T_Element
        external PIB_1D_T_Element
!
        write(*,*)' What is the value of the particle mass M (must be an integer)?'
        read(*,*) m
        if (m .LE. 0) THEN
         write(*,*) 'An error occured -  mass can not be less than or equal to zero.'
            stop
         end if
!
        write(*,*)' What is the value of the length of the box L (must be an integer)?'
        read(*,*) l

          if (l .LE. 0) THEN
          write(*,*) 'An error occured -  box length can not be less than or equal to zero.'
             stop
          end if

!
       write(*,*)' What is the value of the quantum number of the first eigenstate n1 (must be an integer)?'
       read(*,*) n1

       write(*,*)' What is the value of the quantum number of the second eigenstate n2 (must be an integer)?'
       read(*,*) n2
      !
900 format(1X, 'Kinetic energy matrix element', I5, ',',I5,' is',F12.5,'.')
!
!     Do the printing job.
!

      write(*,900) n1,n2,PIB_1D_T_Element(m,l,n1)
 
!
!	If (n1 == n2) then
!		print*,PIB_1D_T_Element(m,l,n1)
!	else
!		print *, 0
!	endif
!
       End Program prgm_02_01
!
!	(m,l,n1) = (a,b,c)
!
	real FUNCTION PIB_1D_T_Element(a,b,c)
        IMPLICIT NONE
        real :: a,b
        integer :: c
         PIB_1D_T_Element = -(c**2)/(a*2*(b**2))
        RETURN
        end FUNCTION PIB_1D_T_Element
!

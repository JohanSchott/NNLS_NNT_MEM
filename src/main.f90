program main
use readinpm
use acm
use saveacm
implicit none

! Program to perform a single analytic continuation using a spectrum method.
! Several spectrum methods are or will become available:
! LS, truncated SVD, Tikhonov, (MaxEnt) and (stochastic sampling).
! Reads input file ac.inp 

call readinp() !read settings in ac.inp, read Matsubara file and exact real axis file, if exists. 
call acs()  ! performs analytic continuation 
call saveac() !save real axis information and how well attained function fits to Matsubara data. 


!deallocate variables
deallocate(a%w)
deallocate(a%f)
deallocate(a%wn)
if(allocated(a%rotation)) deallocate(a%rotation)
deallocate(a%a,a%ra)
deallocate(a%ma,a%me,a%am)
if(exact) then
   deallocate(a%we,a%ae)
endif

write(*,*) "Program finished"

end program

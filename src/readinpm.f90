module readinpm
!This module reads inputfile ac.inp.
!Settings are stored in variables in module settingm. 
!Input Matsubara data is stored in object a, created in module fz
use fz
use settingm
use openm
use input
use leastsquare
use asym
use savem
use gitversion
implicit none

contains

subroutine readinp()
implicit none
character(len=800) :: fm,fam,fe ! Matsubara file, default spectra file, exact real axis file
logical :: ExtDefaultM ! when use default model, if read it from real file
integer :: first !first Matsubara index to fit to
integer :: rmeshtype
real(kind=16) :: qasym
integer :: nasym

real(kind=16),allocatable :: wam(:) ! default dense energy mesh
real(kind=16),allocatable :: tmp(:),tmp2(:) !two temporary variables
real(kind=16) :: ga,dx
integer :: i,ok,N1,N2


open(54,file="ac.inp",iostat=ok)
read(54,'(a)') 
read(54,'(a)') 
read(54,'(a)') proj
read(54,'(a)') 
read(54,'(a)') fm
read(54,'(a)') 
read(54,*) nasym,qasym
read(54,'(a)') 
read(54,*) first
read(54,'(a)') 
read(54,*) a%M
read(54,'(a)') 
read(54,*) a%N
read(54,'(a)') 
read(54,*) a%wmin,a%wmax
read(54,*)
read(54,*) rmeshtype
read(54,*)
read(54,*)
read(54,*)
read(54,*)
read(54,*)
read(54,*)
read(54,*)
read(54,*) ftype
read(54,*)
read(54,*) sumrule
read(54,*)
read(54,*) solver
read(54,*)
read(54,*) ExtDefaultM
read(54,*)
read(54,'(a)') fam
read(54,*)
read(54,*) mtr
read(54,*)
read(54,*) pr
read(54,*)
read(54,*) a%eim
read(54,'(a)') 
read(54,'(a)') 
read(54,*) exact
read(54,'(a)') 
read(54,'(a)') fe
close(54)

!Write parameters read and used from input file  
open(23,file=trim(proj)//"_info.dat")
write(23,'(a,a)') " git hash number: ", git_revision
write(23,'(a)') "Program version: 2.0"
write(23,*)
write(23,*) '#Parameter values read and used from input file:'
write(23,*) 'fm: ',trim(fm)
write(23,*) 'nasym=',nasym,'qasym=',qasym
write(23,*) 'first=',first
write(23,*) 'M=',a%M
write(23,*) 'N=',a%N
write(23,*) 'wmin=',a%wmin,"wmax=",a%wmax
write(23,*) "rmeshtype=",rmeshtype
write(23,*) 'ftype=',ftype
if(.not. sumrule) write(23,*) "sumrule=",sumrule
write(23,*) 'solver=',solver
if(solver==3 .or. solver==4) then
   write(23,*) 'fam:',trim(fam)
endif
write(23,*) 'mtr=',mtr
write(23,*) 'pr=',pr
write(23,*) "eim=",a%eim
write(23,*) 'exact=',exact
if(exact) then
   write(23,*) 'fe:',trim(fe)
endif
write(23,*) " "
write(23,*) "#Output information:"
write(23,*) " "
close(23)

allocate(a%w(a%N))
allocate(a%f(a%N))
allocate(a%wn(a%M))
allocate(a%a(a%N))
allocate(a%ra(a%N))
allocate(a%ma(a%M))
allocate(a%me(a%M))
allocate(a%am(a%N))
allocate(a%ae(a%N))


a%realasym=getrealasym(fm,nasym,qasym) !find estimates for a and b in Re[G(i*w_n)] \approx b-a/w_n^2
a%sweightm=getsweightm(fm,nasym,qasym) !find estimates for s in Im[G(i*w_n)] \approx -s/w_n
call getIndices(fm,first,a%M,a%wn,a%me) !pick selected Matsubara points.  
open(23,file=trim(proj)//"_info.dat",position="append")
if(ftype==1) then
   write(*,*) "Shift realpart on Matsubara axis by estimated value from: Re[G(i*w_n)] \approx -a/w_n^2 + b for big w_n"
   write(23,*) "Shift realpart on Matsubara axis by estimated value from: Re[G(i*w_n)] \approx -a/w_n^2 + b for big w_n"
   a%me=a%me-a%realasym(2) !shift real part
   !a%me=a%me-0.6906919 !shift real part
endif
write(*,'(a,F14.7)') "From input data, estimate b=",a%realasym(2)
write(23,'(a,F14.7)') "From input data, estimate b=",a%realasym(2)
write(*,'(a,F14.7)') "From input data, estimate a=",a%realasym(1)
write(23,'(a,F14.7)') "From input data, estimate a=",a%realasym(1)
if(ftype==3) then
   write(*,*) "Sum rule enforced by estimate value from: Im[G(i*w_n)] \approx -s/w_n for big w_n"
   write(23,*) "Sum rule enforced by estimate value from: Im[G(i*w_n)] \approx -s/w_n for big w_n"
endif
write(*,'(a,F14.7)') "From input data, estimate s=",a%sweightm
write(23,'(a,F14.7)') "From input data, estimate s=",a%sweightm
write(*,*)
write(23,*)

!real axis 
if(mod(a%N,2)==0) then !even number of real axis points
   if(ftype==0 .or. ftype==1 .or. ftype==3) then
      stop "Odd number of real axis energy points is strongly advised for ferminons"
   endif
else  !odd number of real axis points
   if(ftype==2) then
      stop "Even number of real axis energy points is strongly advised for bosons, if Matusbara w_n=0 is included"
   endif
endif
if(rmeshtype==0) then
   write(23,*) "Linear mesh created"
   a%w(:)= [(a%wmin+(a%wmax-a%wmin)*(i-1)/(a%N-1),i=1,a%N) ]
elseif(rmeshtype==1) then
   write(23,*) "Logaritmic mesh created, using space scalig parameter gamma=0.5"
   ga=0.5q0
   if(ftype==4) then
      dx=(log(a%wmax-a%wmin+ga)-log(ga))/(a%N-1)
      do i=1,a%N
         a%w(i)=exp( (i-1)*dx+log(ga))-ga+a%wmin
      enddo
      !a%w(1)=10**(-20q0)
   else
      N2=ceiling(a%N*a%wmax/(a%wmax-a%wmin))
      N1=a%N-N2
      dx=(log(-a%wmin+ga)-log(ga))/N1
      do i=2,N1+1
         a%w(N1+2-i)=-(exp( (i-1)*dx+log(ga))-ga)
      enddo
      dx=(log(a%wmax+ga)-log(ga))/(N2-1)
      do i=1,N2
         a%w(N1+i)=exp( (i-1)*dx+log(ga))-ga
      enddo
   endif
else
   stop "rmeshtype input error.."
endif
if(ftype==2) then
   if(a%w(1)>0 .or. a%w(a%N)<0 ) then
      stop "Both negative and positive frequencies are needed"
   endif
elseif(ftype==4) then
   if(a%w(1)<0) then
      write(*,*) "w(1)=",a%w(1)
      stop "Only positive frequencies should be used for ftype==4" 
   endif
endif
!integration weight
a%f(2:a%N-1)=(a%w(3:a%N)-a%w(1:a%N-2))/2
a%f(1)=(a%w(2)-a%w(1))/2
a%f(a%N)=(a%w(a%N)-a%w(a%N-1))/2

close(23)

if(exact) then
  call openf(fe,a%we,tmp,a%ae)
  a%ae=-1/pi*a%ae
endif
if(solver==3 .or. solver==4) then !Determine Default model on real axis
   if(ExtDefaultM) then
      !interpolate default model to points a%w and save in a%am
      call openf(fam,wam,tmp,tmp2)
      call interp(wam,-1/pi*tmp2,a%w,a%am)
   else
      if(ftype==4) then
         a%am=a%realasym(1)/(a%wmax**2-a%wmin**2)
      else
         a%am=a%sweightm/sum(a%f)
      endif
   endif
endif

end subroutine


end module

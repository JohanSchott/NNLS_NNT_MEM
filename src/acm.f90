module acm
!This module's purpose is to do one analytic continuation.
!It uses object a, defined in module fz of type f, for storing one function's values on real and Matsubara axis.
!It performs the continuation according to the settings in module settingm
use fz
use settingm
use tikhonovm
use maxentm
use savem
implicit none

contains

subroutine acs()
implicit none
call ac(a%w,a%f,a%am,a%wn,a%me,a%realasym(1),a%sweightm,a%a,a%ma,a%errs,a%errm) !returns spectrum rho,corresponding Matsubara value and errors errs and errm
end subroutine

subroutine ac(w,f,am,wn,me,realasym,sweight,rho,ma,errs,errm)
implicit none
real(kind=16),intent(in) :: w(:),f(:),am(:),wn(:)
complex(kind=16),intent(in) :: me(:)
real(kind=16),intent(in) :: realasym,sweight
real(kind=16),allocatable,intent(inout) :: rho(:)
complex(kind=16),allocatable,intent(inout) :: ma(:)
real(kind=16),allocatable,intent(out) :: errs(:)
real(kind=16),intent(out) :: errm

call rac(w,f,am,wn,me,realasym,sweight,rho,errs)
call hilbert(w,f,rho,wn,ma) !calculate G(i*w_n,rho(w)) from hilbert transform
errm=sqrt(sum(abs(ma-me)**2))
end subroutine

subroutine rac(w,f,am,wn,me,realasym,sweight,rho,errs) !returns spectrum rho and errs
implicit none
real(kind=16),intent(in) :: w(:),f(:),am(:),wn(:)
complex(kind=16),intent(in) :: me(:)
real(kind=16),intent(in) :: realasym,sweight
real(kind=16),allocatable,intent(inout) :: rho(:)
real(kind=16),allocatable,intent(out) :: errs(:)

real(kind=16),allocatable :: k(:,:),b(:)
real(kind=16),allocatable :: x(:),xam(:),x0(:) !temporary variables
real(kind=16) :: tolx,tolf,tolg !convergence tolerances for MaxEnt
integer :: i,j,n,m

real(kind=16),allocatable :: kc(:,:),bc(:)

character(len=800) :: str

n=size(w,1)
m=size(wn,1)

!Print info about mathematical problem to solve to this file
open(23,file=trim(proj)//"_info.dat",position="append")
if(ftype==0) then
   write(23,*) "Function is set to be Fermionic"
   if(sumrule) then
      write(23,*) "sumrule: int rho(w)dw=1 is an explicit condition"
   else
      write(23,*) "sumrule: int rho(w)dw=1 is NOT an explicit condition"
   endif
elseif(ftype==1) then
   write(23,*) "Function is set to be Self-energy"
   write(23,*) "Shift of real part is done before using Hilbert transform"
elseif(ftype==2) then
   write(23,*) "Function is set to be Susceptibility"
   write(23,*) "No shift of real part is done before using Hilbert transform"
elseif(ftype==3) then
   write(23,*) "Function is set to be weird Fermion"
   write(23,*) "This means sumrule value to fit to is estimated by looking at the asymptotics."
   write(23,*) "It also means that do an additional fitting, to parameter a in -a/w_n^2,&
   where a is estimated by looking at the asymptotics"
elseif(ftype==4) then
   write(23,*) "Function is set to be symmetric Susceptibility"
   write(23,*) "Only fitting to Re[\chi(i*\omega_n)], using kernel:"
   write(23,*) "Re[G(i*\omega_n)]= \int_{wmin}^{wmax} -2\omega/(\omega^2+\omega_n^2) \rho(\omega) d\omega"
   write(23,*) "where wmin should positive (or zero)."
   write(23,*) "No shift of real part is done before using Hilbert transform"
else
   write(23,*) "Wrong ftype is entered"
   write(23,*) "No shift of real part is done before using Hilbert transform"
endif

!Build hilbert matrix that define analytic continuation problem
if(ftype==0) then
   if(sumrule) then
      allocate(k(2*m+1,n))
   else
      allocate(k(2*m,n))
   endif
elseif(ftype==1 .or. ftype==2) then
   allocate(k(2*m,n))
elseif(ftype==3) then
   allocate(k(2*m+2,n))
elseif(ftype==4) then
   allocate(k(m,n))
endif

if(size(k,1)<size(k,2)) stop "Fewer equations than variables"
write(23,*) "Fit to ",m," Matsubara points"
write(23,*) "size of matrix k, in k*rho=b, is (",size(k,1),",",size(k,2),")" 

if(ftype==4) then
   do i=1,m
      do j=1,n
         k(i,j)=2q0*f(j)*(-w(j))/(wn(i)**2q0+w(j)**2q0)
      enddo
   enddo
else
   do i=1,m
      do j=1,n
         k(i,j)=f(j)*(-w(j))/(wn(i)**2q0+w(j)**2q0)
         k(m+i,j)=f(j)*(-wn(i))/(wn(i)**2q0+w(j)**2q0)
      enddo
   enddo
endif
if(ftype==0 .and. sumrule) then
   k(2*m+1,:)=f    !Im[G(i*w_n)] \approx -s/w_n. fit to s-value
elseif(ftype==3) then
   k(2*m+1,:)=f    !Im[G(i*w_n)] \approx -s/w_n. fit to s-value
   k(2*m+2,:)=f*w  !Re[G(i*w_n)] \approx -a/w_n^2 + b. fit to a-value
endif
!If ftype==2, function is boson, eg. susceptibility, reformulate rho(w)=sign(w)*rho(w), so should have only positive values.
!This is done by transforming k=sign(w)*k and using sign(w)*am(:) as default model.
allocate(xam(n))
xam=am
if(ftype==2) then !boson
   do i=1,n
      if(w(i)<0) then
         k(:,i)=-k(:,i)
         xam(i)=-am(i)
      endif
   enddo
endif
!Matsubara input data and eventually extra constrains
if(ftype==0) then
   if(sumrule) then
      allocate(b(2*m+1))
   else
      allocate(b(2*m))
   endif
elseif(ftype==1 .or. ftype==2) then
   allocate(b(2*m))
elseif(ftype==3) then
   allocate(b(2*m+2))
elseif(ftype==4) then
   allocate(b(m))
endif

if(ftype==4) then
   b(1:m)=real(me) 
else
   b(1:m)=real(me) 
   b(m+1:2*m)=aimag(me)
endif
if(ftype==0 .and. sumrule) then
   b(2*m+1)=1q0
elseif(ftype==3) then
   b(2*m+1)=sweight
   b(2*m+2)=realasym
endif


! Test to add more weight on fitting to lower Matsubara points
!do i=1,5
!   b(i)=10*b(i)
!   b(m+i)=10*b(m+i)
!   k(i,:)=10*k(i,:)
!   k(m+i,:)=10*k(m+i,:)
!enddo
!Test to make constrain int rho(w)dw=1 weaker
!k(2*m+1,:)=0.0001*k(2*m+1,:)
!b(2*m+1)=0.0001*b(2*m+1)


allocate(kc(size(k,1),size(k,2)),bc(size(b,1))) !copies
kc=k
bc=b

if(solver==1) then 
   if(ftype==0 .or. ftype==1 .or. ftype==2 .or. ftype==4) then
      write(23,*) "Solve: min_{x>0} |k*x-b|^2"
      call nnls(k,b,mtr,pr,x) !min_{x>0} |k*x-b|^2  
   elseif(ftype==3) then
      write(23,*) "Solve: min_{x} |k*x-b|^2"
      call ls(k,b,mtr,pr,x) !min_x |k*x-b|^2
   endif
elseif(solver==2) then
   if(ftype==0 .or. ftype==1 .or. ftype==2 .or. ftype==4) then
      write(23,*) "Solve: min_{x>0} |k*x-b|^2 + |alpha*x|^2"
      call tikhonov(k,b,mtr,pr,.true.,x) !min_{x>0} |k*x-b|^2 + |alpha*x|^2 
   elseif(ftype==3) then
      write(23,*) "Solve: min_{x} |k*x-b|^2 + |alpha*x|^2"
      call tikhonov(k,b,mtr,pr,.false.,x) !min_x |k*x-b|^2 + |alpha*x|^2 
   endif
elseif(solver==3) then
   if(ftype==0 .or. ftype==1 .or. ftype==2 .or. ftype==4) then
      write(23,*) "Solve: min_{x>0} |k*x-b|^2 + |alpha*(x-xam)|^2"
      call tikhonov(k,b,xam,mtr,pr,.true.,x) !min_{x>0} |k*x-b|^2 + |alpha*(x-xam)|^2 
   elseif(ftype==3) then
      write(23,*) "Solve: min_{x} |k*x-b|^2 + |alpha*(x-xam)|^2"
      call tikhonov(k,b,xam,mtr,pr,.false.,x) !min_x |k*x-b|^2 + |alpha*(x-xam)|^2 
   endif
elseif(solver==4) then
   allocate(x0(n))
   x0=xam !starting point is equal to default model
   !x0=exp(-w**2)/integ(w,exp(-w**2))*integ(w,xam) !same weight as default model
   tolx=10q0**(-9q0)
   tolf=10q0**(-9q0)
   tolg=10q0**(-9q0)
   if(ftype==0 .or. ftype==1 .or. ftype==2 .or. ftype==4) then
      write(23,*) "Solve: min_{x>0} |k*x-b|^2 + alpha*S[x,xam]"
      call maxent(f,k,b,xam,x0,mtr,pr,.true.,tolx,tolf,tolg,x) !min_{x>0} |k*x-b|^2 + alpha*S[x,xam]  ,starting minimisation with x=x0
   elseif(ftype==3) then
      write(23,*) "Solve: min_{x} |k*x-b|^2 + alpha*S[x,xam]"
      call maxent(f,k,b,xam,x0,mtr,pr,.false.,tolx,tolf,tolg,x) !min_x |k*x-b|^2 + alpha*S[x,xam]     ,starting minimisation with x=x0
   endif
   deallocate(x0)
else
   stop "wrong solver value"
endif

if(allocated(errs)) deallocate(errs)
allocate(errs(3))
errs(1)=sqrt(sum(abs(matmul(k,x)-b)**2)) ! |k*x-b|
deallocate(k,b,xam)
if(allocated(rho)) deallocate(rho)
allocate(rho(n))
rho=x !copy from solver
deallocate(x)
!Transform back if have boson, rho(w)=sign(w)*rho(w)
if(ftype==2) then
   do i=1,n
      if(w(i)<0) then
         rho(i)=-rho(i)
      endif
   enddo
endif
errs(2)=sqrt(sum(abs(rho-am)**2)) ! |rho-am|
errs(3)=sqrt(sum(abs(rho)**2))  !|rho|
close(23)
end subroutine

subroutine hilbert(w,f,rho,wn,g) !calculate G(i*w_n,rho(w)) from hilbert transform. f is integration weight
implicit none
real(kind=16),intent(in) :: w(:),f(:),rho(:),wn(:)
complex(kind=16),allocatable,intent(out) :: g(:)

real(kind=16),allocatable :: gr(:)
complex(kind=16),allocatable :: k(:,:)
real(kind=16),allocatable :: kr(:,:)
integer :: i,j,n,m
n=size(w,1)
m=size(wn,1)
if(ftype==4) then
   allocate(kr(m,n))
   do i=1,m
      do j=1,n
         kr(i,j) = -2q0*f(j)*w(j)/(w(j)**2+wn(i)**2)
      enddo
   enddo
   if(allocated(g)) deallocate(g)
   allocate(g(m),gr(m))
   gr=matmul(kr,rho)
   g=cmplx(gr,0q0,kind=16)
   deallocate(gr,kr)
else
   allocate(k(m,n))
   do i=1,m
      do j=1,n
         k(i,j)=f(j)/(cmplx(0q0,wn(i),kind=16)-w(j))
      enddo
   enddo
   if(allocated(g)) deallocate(g)
   allocate(g(m))
   g=matmul(k,rho)
   deallocate(k)
endif
end subroutine

end module

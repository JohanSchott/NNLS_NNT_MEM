module smearm
implicit none

contains

function smearfapr(x,y,xn,eim)
!Approximatly smears function values y at points x by 
!evaluating spectrum at xn+1i*eim
implicit none
real(kind=16),allocatable :: smearfapr(:)
real(kind=16),intent(in) :: x(:),y(:),xn(:),eim

real(kind=16),parameter :: pi=3.14159265359q0
real(kind=16),allocatable :: f(:),k(:,:)
integer :: n,m,i

n=size(x,1)
m=size(xn,1)
allocate(f(n),k(m,n),smearfapr(m))

f(1)=(x(2)-x(1))/2
f(n)=(x(n)-x(n-1))/2
f(2:n-1)=(x(3:n)-x(1:n-2))/2
do i=1,m
   k(i,:)=f/((xn(i)-x)**2+eim**2)
enddo
smearfapr=eim/pi*matmul(k,y)
deallocate(f,k)
end function

function smearf(x,y,xn,eim,ftype)
!smears function values y at points x by 
!evaluating spectrum at xn+1i*eim using Hilbert transform
!If ftype==4, all x should be positive and y should be odd. 
!Therefor is in that case x and y swapped before Hilber transform is used.
implicit none
real(kind=16),allocatable :: smearf(:)
real(kind=16),intent(in) :: x(:),y(:),xn(:),eim
integer :: ftype

real(kind=8) :: eimd
real(kind=8),allocatable :: xnd(:),smearfd(:)
real(kind=8),allocatable :: tmp(:),w(:),f(:)
real(kind=8),parameter :: pi=3.14159265359d0
real(kind=8) :: c1,c2,d1,d2
integer :: n,m,i,j

n=size(x,1) !nbr of input points
if(ftype==4) then
   allocate(tmp(n))
   tmp=x
   if(tmp(1)==0d0) then !can't take log of zero
      tmp(1)=10d0**(-7)
   endif
   allocate(w(2*n),f(2*n))
   w(1:n)=-tmp(n:1:-1)
   w(n+1:2*n)=tmp
   f(1:n)=-y(n:1:-1)
   f(n+1:2*n)=y
   n=2*n
   deallocate(tmp)
else
   allocate(w(n),f(n))
   w=x
   f=y
endif
m=size(xn,1) !nbr of output points 
allocate(smearf(m),smearfd(m))
allocate(xnd(m))
xnd=xn
eimd=eim

smearfd=0d0
do j=1,m !loop over all output points xn(j)
   do i=1,n-1 !loop over all input intervals: from w(i) to w(i+1) 
      c1=(f(i+1)-f(i))/(w(i+1)-w(i))
      d1=1d0/2*log(((w(i+1)-xnd(j))**2+eimd**2)/((w(i)-xnd(j))**2+eimd**2)) + &
      xnd(j)/eimd*( atan((w(i+1)-xnd(j))/eimd)-atan((w(i)-xnd(j))/eimd) )
      
      c2=(f(i)*w(i+1)-f(i+1)*w(i))/(w(i+1)-w(i))
      d2=( atan((w(i+1)-xnd(j))/eimd)-atan((w(i)-xnd(j))/eimd) )/eimd
      smearfd(j)=smearfd(j)+c1*d1+c2*d2
   enddo
enddo
deallocate(w,f)
deallocate(xnd)
smearf=eimd/pi*smearfd
!write(*,*) "smearf(1)=",smearf(1)
!return smearf
end function

end module

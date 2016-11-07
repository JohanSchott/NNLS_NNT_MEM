module maxentm
use leastsquare
use eigenm
implicit none

contains

subroutine maxent(f,k,b,xm,x0,mtr,pr,nonegx,tolx,tolf,tolg,x)
! If nonegx==.true. then:
!   min_{x>=0} Q(x,alpha)
! If nonegx==.false. then:
!   min_{x} Q(x,alpha)

! Q(x,alpha) = |k*x-b|^2 + alpha*S[x,xam] 
! with S=int x*ln(x/xm) , where S>0 for x/=xam and S=0 for x=xam

! Returns vector x 

! Loop over different alpha and select optimal alpha according to L-curve technique, namely 
! min_alpha |k*x-b|^2*S
! But to find roughly which alpha values to use, start with an iterative algorithm over alpha.
! An lower bound of 10^-11 for alpha is used for numerical reasons. 

! Minimization of Q(x), for a fixed alpha, uses a modified Newton search algorithm.
! Updates are done according to: x=xold+a*dx
! dx is obtained from solving: hessian(x)*dx=-gradient(x) , where hessian and gradient is of function Q.
! Parameter a is obtained by line-search minimization of: Q(xold+a*dx)
! If line-search gives a turning x=xold+a*dx negative for some elemement, adjust those to:
!   x(i)=10^(-10) if nonegx==.true.

! Convergence is reached when |Q(xold)-Q(x)|/Q(x) < tolf  or |grad(xold)-grad(x)|/|grad(x)| < tolg

! Starting minimisation with point x=x0

implicit none
real(kind=16),intent(in) :: f(:),k(:,:),b(:),xm(:),x0(:)
integer,intent(in) :: mtr,pr
logical,intent(in) :: nonegx
real(kind=16),intent(in) :: tolx,tolf,tolg !specifies convergence tolerance  
real(kind=16),allocatable :: x(:)

real(kind=16),allocatable :: lam(:,:),eigen(:)
integer :: n,m,i,ia
integer :: ita,na
real(kind=16) :: alpha,alpha_old
real(kind=16),allocatable :: alphav(:),xt(:,:),qv(:),res(:),st(:)
character(len=800) :: str
na=25 !number of alphas
m=size(k,1)
n=size(k,2)
if(allocated(x)) deallocate(x)
allocate(x(n))
x=x0 !starting point in minimization precedure
open(31,file="output_maxent.dat")
write(31,'(a,E13.6,a,E13.6)') "Q(x0)=",Q(x),"   |grad(x0)|=",sqrt(sum(grad(x)**2))
!Find optimal alpha according to Jarrel:
!alpha_new = |k*x-b|^2/n*1/S(x)*sum_i e_i/(e_i+alpha), where e_i is eigenvalue to matrix lam
!lam_{i,j} = sqrt(x_i) d^2L/(dx_i dx_j) sqrt(x_j) , where L is the likelihood function, in our case: L(x) = |k*x-b|^2
!lam_{i,j}= sqrt(x_i) 2 (k^T*k)_{i,j} sqrt(x_j)
allocate(lam(n,n),eigen(n))
write(*,'(a)') "Start optimization of alpha"
write(31,'(a)') "Start optimization of alpha"
alpha=10q0**(-2q0) !starting guess for alpha
alpha_old=100*alpha
write(31,'(a,E10.3)') "alpha=",alpha
write(*,'(a,E10.3)') "alpha=",alpha
ita=0
do while( 5 < max(alpha_old/alpha,alpha/alpha_old) .and. ita < 10 ) !dalpha=max(alpha_old/alpha,alpha/alpha_old)
   ita=ita+1
   call minQ() !update x
   lam=2*matmul(transpose(k),k)
   do i=1,n
      lam(i,:)=lam(i,:)*sqrt(x)
      lam(:,i)=lam(:,i)*sqrt(x)
   enddo
   eigen=eig(lam) !get eigenvalues to lam
   alpha_old=alpha
   !alpha=1q0/(2*entropy(x))*sum(eigen/(eigen+alpha))  !update alpha
   alpha=residue2(x)/m*1q0/(entropy(x))*sum(eigen/(eigen+alpha))  !update alpha
   if(alpha<0q0) then
      alpha=10**(-10q0)
      write(31,'(a)') 'Negative alpha proposed. Use a small alpha instead.'
   endif
   !mixing old and new alpha
   alpha = 10q0**(log(alpha*alpha_old)/(2*log(10q0)))
   if(alpha<1.1*10**(-10q0) ) then
      write(31,'(a)') "alpha smaller than 10^(-10) proposed in the iterative algorithm."
      write(31,'(a)') "alpha put to 10^(-10) instead."
      alpha=10q0**(-10q0)
      exit
   else
      write(31,*) "alpha=",alpha
      write(*,*) "alpha=",alpha
   endif
enddo
deallocate(lam,eigen)

write(31,'(a)') 
write(*,'(a)') "Start exploring L-curve for optimal alpha"
write(31,'(a)') "Start exploring L-curve for optimal alpha"
write(31,'(a,E10.3)') "Search around converged alpha value, which is alpha=",alpha
allocate(alphav(na),xt(n,na),qv(na),res(na),st(na))
st = -1
!alphav= [ (10q0**(-(8q0+(i-1q0)*(12q0-8q0)/(na-1q0))),i=1,na) ]
alphav= [ (alpha*10q0**(-(-3q0+(i-1q0)*(3q0+2q0)/(na-1q0))),i=1,na) ]
ia=0
do ita=1,na !Loop A
   ia=ia+1
   alpha=alphav(ita)
   write(31,'(a,E10.3)') "alpha=",alpha
   write(*,'(a,E10.3)') "alpha=",alpha
   call minQ() !update x
   !save solution for particular alpha
   xt(:,ita)=x
   qv(ita)=Q(x) 
   res(ita)=sqrt(residue2(x)) ! |A*x-b|
   st(ita)=entropy(x)  ! S , entropy
   if (isnan(qv(ita))) return  ! check for NaN 
enddo !end loop A
!count number of alphas giving positive entropy. Want at least two.
i=0
do ita=1,ia
   if(st(ita)>0) i=i+1
enddo

do while(i<2) !If less than two solutions with positive entropy, try smaller alphas
   write(31,'(a)') "Only",i,"solutions with positive entropy found." 
   write(31,'(a)') "Try to use smaller alphas. use alpha from",alpha*10**(1q0)," down to",alpha*10**(-0.5q0)
   alphav= [ (alpha*10q0**(-(-1q0+(i-1q0)*(0.5q0+1q0)/(na-1q0))),i=1,na) ]
   ia=0
   st = -1
   do ita=1,na !Loop A
      ia=ia+1
      alpha=alphav(ita)
      write(31,'(a,E10.3)') "alpha=",alpha
      write(*,'(a,E10.3)') "alpha=",alpha
      call minQ() !update x
      !save solution for particular alpha
      xt(:,ita)=x
      qv(ita)=Q(x) 
      res(ita)=sqrt(residue2(x)) ! |A*x-b|
      st(ita)=entropy(x)  ! S , entropy
      if (isnan(qv(ita))) return  ! check for NaN 
   enddo !end loop A
   !count number of alphas giving positive entropy. Want at least two.
   i=0
   do ita=1,ia
      if(st(ita)>0) i=i+1
   enddo
enddo
write(31,'(I4,a)') i," solutions with positive entropy found."
write(31,'(a)') "Sweep over different alphas is finished. Summary of the sweep:"
write(31,'(a)') "       alpha          Q(x)         |A*x-b|^2        S         |A*x-b|^2*S"
do i=1,na
   write(31,'(5E15.5)') alphav(i),qv(i),res(i)**2,st(i),res(i)**2*st(i)
enddo
str=int2str(na)
str="("//trim(str)//"E13.4)"
open(54,file="output_maxent_a.dat")
do i=1,n
   write(54,trim(str)) xt(i,:)
enddo
close(54)
i=minloc(log(res**2*st),1)
x=xt(:,i)
alpha=alphav(i)
write(31,'(a,E15.5,a,I4)') "Pick alpha=",alpha," index=",i

write(31,*)
write(31,'(a,E15.6,a,E15.6)') "Q(x)=",Q(x),"  |grad(x)|=",sqrt(sum(grad(x)**2))
write(31,'(a,E15.6,a,E15.6)') "Q(xm)=",Q(xm),"  |grad(xm)|=",sqrt(sum(grad(xm)**2))
close(31)
deallocate(alphav,xt,qv,res,st)

contains
   
   subroutine minQ()
   !Finds minimun to Q(x), for a particular alpha, and update x to that point
   implicit none
   integer :: i,j !local variables
   integer :: ib,ic,ibmax,icmax
   real(kind=16) :: a,amax
   real(kind=16),allocatable :: da(:),dx(:),xold(:)
   ibmax=3*n !maximum iteration for loop B
   icmax=200 !maximum iteration for loop C
   allocate(da(1),dx(n),xold(n))
   ib=0
   do i=1,ibmax !Loop B
      xold=x !save x before changing its value
      ib=ib+1
      call ls(hess(x),-grad(x),mtr,pr,dx)
      amax=minval(-x/dx,1,dx<0)
      amax=min(amax,1q0)

      !Method 1: Preform a line search along points x+a*dx, wrt "a"
      !Starting with a=0
      a=0q0
      ic=0
      !write(31,*) "Line-search to find good a-value"
      do j=1,icmax !Loop C
         ic=ic+1
         call ls(DgradDa(x,dx,a),-grad(x+a*dx),mtr,pr,da)
         if(a+da(1)>0 .and. a+da(1)<amax) then
            a=a+da(1)
            !write(31,'(a,E10.4)') "a=",a
            if(da(1)/a< 10q0**(-9q0) ) then
               !write(*,'(a)') "da/a<10^(-10). Line search converged"
               exit
            endif
         elseif(a+da(1)>=amax) then
            a=a+da(1)
            !we will with this new a get negative x+a*dx values, but remove them
            !below
            exit
         else
            write(31,*)
            write(31,*) "Don't update a since a+da<0 which is the wrong direction"
            write(31,*) "ic=",ic
            write(31,*) "a=",a
            write(31,*) "da=",da(1)
            write(31,*) "a+da=",a+da(1)
            write(31,*)
            exit
            !stop "negative da should never happen..."
         endif
      enddo
      write(31,'(I5,a)') ic," iterations in Loop C."

      !Method 4: Use factor times amax in the Newton step.
      !a=0.1*amax

      !Allowed values: 0<a<amax to garantee positive new x
      !Have now calculated an a-value, not neccesarry fulfilling that, so adjust negative x-values to small positive values
      !write(31,'(a,E10.3)') "a to be used:",a

      if(nonegx) then
         where(x+a*dx<=0) 
            x=10q0**(-10q0)
         elsewhere
            x=x+a*dx
         endwhere
      else
         x=x+a*dx
      endif
      !Print convergence information:
      !write(31,'(a,E15.6,a,E15.6)') "Q(x)=",Q(x),"  |grad(x)|=",sqrt(sum(grad(x)**2))
      !write(31,*)
      !write(31,'(a,E15.6,a,E15.6)') "|xold-x|=",sqrt(sum((xold-x)**2))," |xold-x|/|x|=",sqrt(sum((xold-x)**2))/sqrt(sum(x**2))
      !write(31,'(a,E15.6,a,E15.6)') "Q(xold)-Q(x)=",Q(xold)-Q(x),"  (Q(xold)-Q(x))/Q(x)=",(Q(xold)-Q(x))/Q(x)
      !write(*,'(a,E15.6,a,E15.6)') "Q(xold)-Q(x)=",Q(xold)-Q(x),"  (Q(xold)-Q(x))/Q(x)=",(Q(xold)-Q(x))/Q(x)
      !write(31,'(a,E15.6,a,E15.6)') "|grad(xold)-grad(x)|=",sqrt(sum((grad(xold)-grad(x))**2)), &
      !"  |grad(xold)-grad(x)|/|grad(x)|=",sqrt(sum((grad(xold)-grad(x))**2))/sqrt(sum(grad(x)**2))
      !Check if to exit loop B
      if(sqrt(sum((xold-x)**2))/sqrt(sum(x**2)) < tolx ) then !relative change in x 
         write(31,'(a)') "Relative change in x considered converged"
         exit
      elseif(abs(Q(xold)-Q(x))/Q(x) < tolf ) then !relative change in function
         write(31,'(a)') "Relative change in function considered converged"
         exit
      elseif(sqrt(sum((grad(xold)-grad(x))**2))/sqrt(sum(grad(x)**2)) < tolg ) then !relative change in gradient
         write(31,'(a)') "Relative change in gradient considered converged"
         exit
      endif
   enddo !end loop B
   
   write(31,'(I5,a)') ib," iterations in Loop B."
   if(ib==ibmax) write(31,'(a)') "Warning: Minimazation using Newton calculation did not converge"
   write(31,*)
   deallocate(dx,da,xold)
   end subroutine

   function residue2(r)
   !Function returning residue squared: |k*r-b|^2
   implicit none
   real(kind=16),intent(in) :: r(:)
   real(kind=16) :: residue2
   residue2=sum((matmul(k,r)-b)**2)
   end function

   function entropy(r)
   !Function returning the entropy: S[r,xm]= sum_i f_i*r_i*ln(r_i/xm_i)
   implicit none
   real(kind=16),intent(in) :: r(:)
   real(kind=16) :: entropy
   entropy=sum(f*r*log(r/xm))
   end function

   function Q(r)
   !Function returning value of function: |k*r-b|^2+alpha*S[r,xm]
   implicit none
   real(kind=16),intent(in) :: r(:)
   real(kind=16) :: Q
   Q=residue2(r)+alpha*entropy(r)
   end function

   function grad(r)
   !Function returning gradient of function Q defined above, namely returning:
   !grad_i(r)=dQ/dr_i(r)=(2k^T*(k*r-b))_i+alpha*f_i*(ln(r_i/xm_i)+1)
   implicit none
   real(kind=16) :: grad(n)
   real(kind=16),intent(in) :: r(:)
   grad=2*matmul(transpose(k),matmul(k,r)-b)+alpha*f*(log(r/xm)+1)
   end function
   
   function DgradDa(r,dr,c)
   !Function returning derivative of grad(r+c*dr) wrt c, namely returning:
   !DgradDa_i(r,dr,c)=(2k^T*(k*dr))_i+alpha*f_i*dr_i/(r_i+c*dr_i)
   implicit none
   real(kind=16) :: DgradDa(n,1)
   real(kind=16),intent(in) :: r(:),dr(:),c
   DgradDa(:,1)=2*matmul(transpose(k),matmul(k,dr))+alpha*f*dr/(r+c*dr)
   end function
   
   function hess(r)
   !function returning hessian of function Q defined above, namely returning:
   !hess_{i,j}=d^2Q/dr^2_{i,j}=(2k^T*k)_{i,j}+alpha*f_i/r_i*delta_{i,j}, where delta_{i,j}=0 if i/=j and delta_{i,i}=1
   implicit none
   real(kind=16) :: hess(n,n)
   real(kind=16),intent(in) :: r(:)
   integer :: i
   hess=2*matmul(transpose(k),k)
   do i=1,size(k,2)
      hess(i,i)=hess(i,i)+alpha*f(i)/r(i)
   enddo
   end function


end subroutine

end module

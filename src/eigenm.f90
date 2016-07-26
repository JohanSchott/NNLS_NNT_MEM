module eigenm
implicit none

interface eig
   module procedure deig,qeig
end interface eig

contains 

function qeig(A)
!calculate eigenvalues of real symmetric matrix A and return eigenvalues in vector x
implicit none
real(kind=16),intent(in) :: A(:,:)
real(kind=16),allocatable :: qeig(:)

integer :: n
real(kind=8),allocatable :: dA(:,:),e(:)

n=size(A,1)
allocate(dA(n,n),e(n))
dA=A !convert to double precision
e=deig(dA)
if(allocated(qeig)) deallocate(qeig)
allocate(qeig(n))
qeig=e
deallocate(dA,e)
end function

function deig(A)
!calculate eigenvalues of real symmetric matrix A and return eigenvalues in vector x
implicit none
real(kind=8),intent(in) :: A(:,:)
real(kind=8),allocatable :: deig(:)

real(kind=8),allocatable :: k(:,:),tmp(:) !computing variables

!dsytrd variables
character(len=1) :: uplo
integer :: n,lda,lwork,info
real(kind=8),allocatable :: d(:),e(:),tau(:),work(:)

!dstegr variables (in top of those already specified for dsytrd)
character(len=1) :: jobz,range
real(kind=8) :: vl,vu,abstol
integer :: il,iu,m,ldz,liwork
real(kind=8),allocatable :: w(:),z(:,:)
integer,allocatable :: isuppz(:),iwork(:)

n=size(A,1)
allocate(k(n,n))
k=A !copy

!check if matrix is symmetric
if( sum(abs(k-transpose(k)))/(n*n-n) > 10**(-10q0) ) stop "matrix not symmetric, which is needed for dsytrd" 

!factorize matrix k: k=Q*T*transpose(Q), where T is a tridiagonal matrix
!Using lapack routine: dsytrd
uplo="U"
lda=n
allocate(d(n),e(n-1),tau(n-1))
allocate(work(1))
lwork=-1
call dsytrd(uplo,n,k,lda,d,e,tau,work,lwork,info) !get optimal lwork
if(info<0) then
   write(*,*) "Argument",-info," in dsytrd had illegal value"
   stop
endif
lwork=work(1)
deallocate(work)
allocate(work(lwork))
call dsytrd(uplo,n,k,lda,d,e,tau,work,lwork,info)
deallocate(work)
!tridiagonal elements are stored in arrays: d (diagonal) and e (subdiagonal)

!find eigenvalues of tridiagonal matrix T 
!Using lapack routine: dstegr
!Pass diagonal and subdiagonal of T in d and e
jobz="N"
range="A"
allocate(tmp(n-1))
tmp=e
deallocate(e)
allocate(e(n))
e(1:n-1)=tmp
deallocate(tmp)
allocate(w(n))
ldz=1
allocate(z(ldz,n))
allocate(isuppz(n))
allocate(work(1),iwork(1))
lwork=-1
liwork=-1
call dstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
if(info<0) then
   write(*,*) "Argument",-info," in dstegr had illegal value"   
   stop
elseif(info>0) then
   if(10<=info .and. info<20) then
      stop "dstegr stopped due to internal error in DLARRE"
   elseif(20<=info .and. info<30) then
      stop "dstegr stopped due to internal error in DLARRV"
   else
      stop "dstegr stopped unexpectedly"
   endif
endif
lwork=work(1)
liwork=iwork(1)
deallocate(work,iwork)
allocate(work(lwork),iwork(liwork))
call dstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
if(m/=n) stop "Not all eigenvalues found in dstegr"
deallocate(d,e,work,iwork,z,isuppz)
!eigenvalues are calculated and returned to variable w

if(allocated(deig)) deallocate(deig) !prepare return array
allocate(deig(n))
deig=w !returns eigvalues
deallocate(w)
deallocate(k)
end function

end module

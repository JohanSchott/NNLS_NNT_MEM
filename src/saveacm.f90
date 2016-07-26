module saveacm
!This module's purpose is to save everyting in the object a, defined in module fz of type f.
use fz
use savem
use settingm
use matlab
use smearm
implicit none


contains

subroutine saveac()
implicit none
character(len=800) :: f
real(kind=16),allocatable :: tmp(:),tmp2(:),tmp3(:) !temporary variables

!save G(i*w_n)[rho]
f=trim(proj)//"_mats.dat"
call save2f(f,a%wn,real(a%ma),aimag(a%ma),real(a%me),aimag(a%me))

!save spectrum rho(w)
f=trim(proj)//"_a.dat"
if(solver==3 .or. solver==4) then
   call save2f(f,a%w,a%a,a%am)
else
   call save2f(f,a%w,a%a)
endif
!save smeared rho(w)
f=trim(proj)//"_smear.dat"
if(exact) then
   allocate(a%as(size(a%ae,1)))
   a%as=smearf(a%w,a%a,a%we,a%eim,ftype)
   call save2f(f,a%we,a%as)
else
   allocate(a%as(a%n))
   a%as=smearf(a%w,a%a,a%w,a%eim,ftype)
   call save2f(f,a%w,a%as)
endif
!Print devation values between smeared and exact spectrum 
if(exact) then
   f=trim(proj)//"_info.dat"
   open(23,file=trim(f),position="append")
   write(23,*)
   write(23,'(a)') "Print devation values between smeared and exact spectrum:" 
   write(23,'(a,E13.4)') "int |Asmear(w)-Aexact(w)|dw=",integ(a%we,abs(a%as-a%ae))
   write(23,'(a,E13.4)') "int |Asmear(w)-Aexact(w)|dw / int |Aexact(w)|dw=",integ(a%we,abs(a%as-a%ae))/integ(a%we,abs(a%ae))
   write(23,'(a,E13.4)') "int (Asmear(w)-Aexact(w))^2dw=",integ(a%we,(a%as-a%ae)**2)
   write(23,'(a,E13.4)') "int (Asmear(w)-Aexact(w))^2dw / int (Aexact(w))^2dw=",integ(a%we,(a%as-a%ae)**2)/integ(a%we,(a%ae)**2)
   write(23,*)
   close(23)
endif
deallocate(a%as)


!Print deviation values
f=trim(proj)//"_info.dat"
open(23,file=trim(f),position="append")
write(23,'(a,E13.4)') "Solver error: |k*x-b|=",a%errs(1)
write(*,'(a,E13.4)') "Solver error: |k*x-b|=",a%errs(1)
if(ftype==0) then
   write(23,'(a,E13.4)') "Check, should be same as above: &
   sqrt(|G[rho]-G|^2+|int rho(w)dw-1|^2)=",sqrt(a%errm**2+(integ(a%w,a%a)-1)**2)
   write(*,'(a,E13.4)') "Check, should be same as above:  &
   sqrt(|G[rho]-G|^2+|int rho(w)dw-1|^2)=",sqrt(a%errm**2+(integ(a%w,a%a)-1)**2)
elseif(ftype==1 .or. ftype==2) then
   write(23,'(a,E13.4)') "Check, should be same as above: |G[rho]-G|=",a%errm
   write(*,'(a,E13.4)') "Check, should be same as above: |G[rho]-G|=",a%errm
elseif(ftype==3) then
   write(23,'(a,E13.4)') "Check, should be same as above: &
   sqrt(|G[rho]-G|^2+|int rho(w)dw-s|^2+|int w*rho(w)-a|^2)=",a%errm
   write(*,'(a,E13.4)') "Check, should be same as above:  &
   sqrt(|G[rho]-G|^2+|int rho(w)dw-s|^2+|int w*rho(w)-a|^2)=",a%errm
elseif(ftype==4) then
   write(23,'(a,E13.4)') "Check, should be same as above: |Re[G[rho]-G]]|=",a%errm
   write(*,'(a,E13.4)') "Check, should be same as above: |Re[G[rho]-G]|=",a%errm
endif
write(23,'(a,E13.4)') "Matsubara error: |G[rho]-G)|=",a%errm
write(*,'(a,E13.4)') "Matsubara error: |G[rho]-G)|=",a%errm
!Print info about asymptotics and spectral weight
if(ftype==0) then
   write(23,'(a,E13.4)') "|int rho(w)dw-1|=",abs(integ(a%w,a%a)-1)
   write(*,'(a,E13.4)') "|int rho(w)dw-1|=",abs(integ(a%w,a%a)-1)
elseif(ftype==3) then
   write(23,'(a,E13.4)') "|int rho(w)dw-s|=",abs(integ(a%w,a%a)-a%sweightm)
   write(*,'(a,E13.4)') "|int rho(w)dw-s|=",abs(integ(a%w,a%a)-a%sweightm)
   write(23,'(a,E13.4)') "|int w*rho(w)dw-a|=",abs(integ(a%w,a%w*a%a)-a%realasym(1))
   write(*,'(a,E13.4)') "|int w*rho(w)dw-a|=",abs(integ(a%w,a%w*a%a)-a%realasym(1))
endif
if(solver==3 .or. solver==4) then
   write(23,'(a,E13.4)') "|rho-am|=",a%errs(2)
   write(*,'(a,E13.4)') "|rho-am|=",a%errs(2)
endif
write(23,'(a,E13.4)') "|rho|=",a%errs(3)
write(*,'(a,E13.4)') "|rho|=",a%errs(3)
write(23,'(a,E13.4)') "int rho(w) dw=",integ(a%w,a%a)
write(*,'(a,E13.4)') "int rho(w) dw=",integ(a%w,a%a)
write(*,*)
close(23)


!save solver specific output
if(solver==2 .or. solver==3) then !Tikhonov
   f=trim(proj)//"_tikhonov.dat"
   call system("mv output_tikhonov.dat "//trim(f))
   f="tmp.dat"
   call save2f(f,a%w)
   f=trim(proj)//"_tikhonov_a.dat"
   call system("paste tmp.dat output_tikhonov_a.dat > "//trim(f))
   call system("rm tmp.dat output_tikhonov_a.dat")
elseif(solver==4) then !MaxEnt
   f=trim(proj)//"_maxent.dat"
   call system("mv output_maxent.dat "//trim(f))
   f="tmp.dat"
   call save2f(f,a%w)
   f=trim(proj)//"_maxent_a.dat"
   call system("paste tmp.dat output_maxent_a.dat > "//trim(f))
   call system("rm tmp.dat output_maxent_a.dat")
endif

end subroutine

end module 

module settingm

implicit none

character(len=800) :: proj ! Simulation name
logical :: exact           ! If real axis solution exists
integer :: ftype           ! 0:Fermionic Green's function, 1:Self-energy, 2:Susceptibility, 3:Weird Fermion and 4:odd spectral susceptibility. 
logical :: sumrule         ! Only used for ftype==0. If .false., s=1 is not a condition.
integer :: solver          ! 1 for LS, 2 for Tikhonov with min |rho|^2, 3 for Tikhonov with min |rho-am|^2, 4 for MEM with entropy S=int rho*ln(rho/am)dw
integer :: mtr             ! If (0) 1, (no) multiplication with transpose matrix. 2 means explicit SVD. 3 means truncated SVD.
integer :: pr              ! Which inversion precision to use. 64: double, 128: quadruple
integer :: covariance     ! 0:no covariance is read. 1:the cholesky matrix W from covariance^-1 = W^T*W is read.

real(kind=16),parameter :: pi=3.14159265359

contains

end module

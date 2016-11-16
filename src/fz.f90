module fz

implicit none

type f ! Describes a function in the complex plane.
   integer                      :: N             ! number of spectral points.
   integer                      :: M             ! number of matsubara points to fit to.
   real(kind=16)                :: wmin          ! minimum value of real-axis mesh
   real(kind=16)                :: wmax          ! maximum value of real-axis mesh
   real(kind=16),allocatable    :: w(:)          ! real-axis energies
   real(kind=16),allocatable    :: f(:)          ! integration weights
   real(kind=16),allocatable    :: wn(:)         ! input points
   real(kind=16)                :: eim           ! smearing on real axis
   real(kind=16),allocatable    :: a(:)          ! spectral function 
   real(kind=16),allocatable    :: as(:)         ! smeared spectral function
   complex(kind=16),allocatable :: ra(:)         ! real axis
   complex(kind=16),allocatable :: ma(:)         ! using the spectral function, the values obtained at the input points
   complex(kind=16),allocatable :: me(:)         ! input values
   real(kind=16),allocatable    :: am(:)         ! default spectral function
   real(kind=16),allocatable    :: ae(:)         ! exact spectral function
   real(kind=16),allocatable    :: we(:)         ! energy associated with read exact spectral function
   real(kind=16),allocatable    :: errs(:)       ! errors in solving matrix problem
   real(kind=16)                :: errm,errx     ! error on matsubara axis and on real axis
   real(kind=16)                :: sweight       ! spectral weight. calculated from: int A(w) dw 
   real(kind=16)                :: sweightm      ! spectral weight. calculated from: Im[G(i*w_n)] \approx -s/w_n
   real(kind=16)                :: realasym(2)   ! Re[G(i*w_n)] \approx  -realasym(1)/w_n^2 + realasym(2)
   real(kind=16),allocatable    :: rotation(:,:) ! the rotation matrix from the covariance matrix. Has to have dimensions (M,M)
end type f

type(f) :: a

contains

end module

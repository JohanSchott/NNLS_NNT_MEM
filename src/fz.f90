module fz

implicit none

type f ! Describes a function in the complex plane.
   integer :: N,M  ! number of spectral points and number of matsubara points to fit to
   real(kind=16) :: wmin,wmax !specify real axis window
   real(kind=16),allocatable :: w(:),f(:) !energies and integration weights
   real(kind=16),allocatable :: wn(:) !Matsubara points
   real(kind=16) :: eim  !smearing on real axis
   real(kind=16),allocatable :: a(:),as(:) !spectral function and smeared spectral function
   complex(kind=16),allocatable :: ra(:),ma(:),me(:) ! real axis, Matsubara axis values and exact Matsubara.
   real(kind=16),allocatable :: am(:) ! default spectral function
   real(kind=16),allocatable :: we(:),ae(:) ! energy and exact spectral function
   real(kind=16),allocatable :: errs(:) ! errors in solving matrix problem
   real(kind=16) :: errm,errx ! error on matsubara axis and on real axis
   real(kind=16) :: sweight,sweightm !spectral weight. calculated from: int A(w) dw and by Im[G(i*w_n)] \approx -s/w_n
   real(kind=16) :: realasym(2) ! Re[G(i*w_n)] \approx  -realasym(1)/w_n^2 + realasym(2)
end type f


type(f) :: a

contains

end module

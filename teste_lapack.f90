SUBROUTINE MKL_test
USE lapack95

implicit none

real aaa(3,3),bbb(3)

integer piv(3)

data aaa/2.5,-1.0,-1.0, &

      -1.0, 2.5,-1.0, &

      -1.0,-1.0, 2.5/

data bbb/1.0, 2.0, 3.0/

!

call getrf(aaa,piv)

!

call getrs(aaa,piv,bbb)

write(*,*)bbb



end SUBROUTINE MKL_test
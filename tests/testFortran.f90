PROGRAM Lorenz63_NILSS

    IMPLICIT NONE

    REAL(8) :: x(2), y(3, 2)
    REAL(8), ALLOCATABLE :: z(:,:)

    INTEGER :: i

    x(1) = 1.0
    x(2) = 2.0

    DO i = 1, 3
        y(i, 1) = i
        y(i, 2) = i + 10
    END DO

    ALLOCATE(z(3,2))
    z(:,:) = y(:,:)

    CALL c_subroutine(x, y, z)

END PROGRAM

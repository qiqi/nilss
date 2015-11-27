PROGRAM Lorenz63_NILSS

    IMPLICIT NONE

    INTEGER, PARAMETER :: nSteps = 10000000

    REAL(8) :: s(1)
    REAL(8) :: x(3)
    REAL(8) :: z
    INTEGER :: iStep

    s(1) = 25.0
    WRITE (*,'(A4,X,A16)') 'S', 'Z'
    DO WHILE (s(1) < 40)
        x(:) = 1.0
        x(2) = 10.0
        DO iStep = 1, 10000
            CALL Step(x, s)
        END DO

        DO iStep = 1, nSteps
            CALL Step(x, s)
            z = z + x(3)
        END DO
        z = z / nSteps

        WRITE (*,'(F4.1,X,F16.12)') s, z

        s(1) = s(1) + 1
    END DO

CONTAINS

SUBROUTINE Step(x, s)
    IMPLICIT NONE
    REAL(8), INTENT(inout) :: x(3)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.001

    REAL(8) :: dx(3)
    dx(1) = 10 * (x(2) - x(1))
    dx(2) = x(1) * (s(1) - x(3)) - x(2)
    dx(3) = x(1) * x(2) - 8. / 3 * x(3)
    x(:) = x(:) + dt * dx(:)
END SUBROUTINE

END PROGRAM

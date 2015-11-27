PROGRAM Lorenz95_NILSS

    IMPLICIT NONE

    INTEGER, PARAMETER :: nState = 40
    INTEGER, PARAMETER :: nSteps = 10000000

    REAL(8) :: s(1)
    REAL(8) :: x(nState)
    REAL(8) :: xSum
    INTEGER :: iStep

    s(1) = 20.0
    WRITE (*,'(A4,X,A16)') 'S', 'X_SUM'
    DO WHILE (s(1) < 40)
        x(:) = 0.0
        x(1) = 10.0
        DO iStep = 1, 10000
            CALL Step(x, s)
        END DO

        DO iStep = 1, nSteps
            CALL Step(x, s)
            xSum = xSum + SUM(x)
        END DO
        xSum = xSum / nSteps

        WRITE (*,'(F4.1,X,F16.12)') s, xSum

        s(1) = s(1) + 1
    END DO

CONTAINS

SUBROUTINE Step(x, s)
    IMPLICIT NONE
    REAL(8), INTENT(inout) :: x(nState)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.001
    INTEGER :: i, ip, im, imm

    REAL(8) :: dx(nState)
    DO i = 1, nState
        ip = mod(i, nState) + 1
        im = mod(i + nState - 2, nState) + 1
        imm = mod(i + nState - 3, nState) + 1
        dx(i) = s(1) - x(imm) * x(im) + x(im) * x(ip) - x(i)
    END DO
    x(:) = x(:) + dt * dx(:)
END SUBROUTINE

END PROGRAM

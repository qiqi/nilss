PROGRAM Lorenz63_NILSS

    IMPLICIT NONE

    REAL(8), PARAMETER :: s(1) = (/30.0/)
    INTEGER, PARAMETER :: nChunks = 500, nStepsPerChunk = 100
    INTEGER, PARAMETER :: nHomoAdjoint = 2

    INTEGER :: iStep, iChunk, iAdjoint
    REAL(8) :: x(3), y(3, nHomoAdjoint+1), dotProductWeights(3)
    REAL(8) :: grad(1, nHomoAdjoint+1), lss_grad(1)
    REAL(8), ALLOCATABLE :: history(:,:,:)

    INTEGER :: nTotalSteps
    REAL(8) :: windowFunction, timeFraction

    WRITE (*,*) 'TESTING ADJOINT IMPLEMENTATION'
    DO iChunk = 1, 10
        CALL TestAdjoint
    END DO

    WRITE (*,*) 'INITIAL TRANSIENT'
    x(:) = (/1.0, 1.0, 1.0/)
    DO iStep = 1, 1000
        CALL Step(x, s)
    END DO

    WRITE (*,*) 'PRIMAL CALCULATION'
    Allocate(history(3, nStepsPerChunk, nChunks))
    DO iChunk = 1, nChunks
        DO iStep = 1, nStepsPerChunk
            history(:, iStep, iChunk) = x(:)
            CALL Step(x, s)
        END DO
    END DO

    WRITE (*,*) 'ADJOINT CALCULATION'
    dotProductWeights(:) = (/1.0, 1.0, 1.0/)
    ! NiLSS_init(nHomoAdjoint, nStateVariables,
    !            nDesignVariables, dotProductWeights)
    CALL NiLSS_init(nHomoAdjoint, 3, 1, dotProductWeights)

    y(:,:) = 0.0
    y(1,1) = 1.0 ! make the terminal y(i,:) different for every i
    y(2,2) = 1.0

    DO iChunk = nChunks, 1, -1
        DO iStep = nStepsPerChunk, 1, -1
            x(:) = history(:, iStep, iChunk)
            DO iAdjoint = 1, nHomoAdjoint + 1
                CALL gradContribution(x, y(:,iAdjoint), s, grad(:,iAdjoint))
                CALL Adjoint(x, y(:,iAdjoint), s)
            END DO
            ! windowing
            nTotalSteps = nChunks * nStepsPerChunk
            timeFraction = real(iStep + iChunk * nStepsPerChunk) / nTotalSteps
            windowFunction = sin(timeFraction) * sin(timeFraction) * 2
            y(3, nHomoAdjoint + 1) = y(3, nHomoAdjoint + 1) + &
                windowFunction / nTotalSteps
            ! checkpoint the adjoint solution
            CALL NiLSS_checkpoint(y, grad)
        END DO
    END DO

    CALL NiLSS_gradient(lss_grad)
    Write (*,*) 'NI-LSS Gradient = ', lss_grad(1)

CONTAINS

SUBROUTINE Step(x, s)
    IMPLICIT NONE
    REAL(8), INTENT(inout) :: x(3)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.005

    REAL(8) :: dx(3)
    dx(1) = 10 * (x(2) - x(1))
    dx(2) = x(1) * (s(1) - x(3)) - x(2)
    dx(3) = x(1) * x(2) - 8. / 3 * x(3)
    x(:) = x(:) + dt * dx(:)
END SUBROUTINE

SUBROUTINE Adjoint(x, y, s)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x(3)
    REAL(8), INTENT(inout) :: y(3)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.005

    REAL(8) :: dy(3)
    dy(1) = -10 * y(1) + (s(1) - x(3)) * y(2) + x(2) * y(3)
    dy(2) = 10 * y(1) - y(2) + x(1) * y(3)
    dy(3) = -x(1) * y(2) - 8./3 * y(3)
    y(:) = y(:) + dt * dy(:)
END SUBROUTINE

SUBROUTINE gradContribution(x, y, s, grad)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x(3), y(3)
    REAL(8), INTENT(in) :: s(1)
    REAL(8), INTENT(inout) :: grad(1)

    REAL(8), PARAMETER :: dt = 0.005

    grad(1) = grad(1) + x(1) * y(2) * dt
END SUBROUTINE

SUBROUTINE testAdjoint
    IMPLICIT NONE
    REAL(8) :: x(1000,3), xp(3), xm(3), y(3), s(1), sp(1), sm(1), eps, grad1, grad2(1)
    INTEGER :: iStep

    CALL RANDOM_NUMBER(x(1,:))
    CALL RANDOM_NUMBER(y)
    CALL RANDOM_NUMBER(s)
    x(1,:) = x(1,:) * 50
    s(:) = s(:) * 100

    eps = 1E-9
    sp(:) = s(:) + eps / 2
    sm(:) = s(:) - eps / 2

    xp(:) = x(1,:)
    xm(:) = x(1,:)
    DO iStep = 2,100
        x(iStep,:) = x(iStep-1,:)
        CALL Step(x(iStep,:), s)
        CALL Step(xp, sp)
        CALL Step(xm, sm)
    END DO
    grad1 = SUM((xp - xm) * y) / eps

    grad2(1) = 0
    DO iStep = 99,1,-1
        CALL gradContribution(x(iStep,:), y, s, grad2)
        CALL Adjoint(x(iStep,:), y, s)
    END DO

    write (*,*) "FD = ", grad1, " ADJ = ", grad2, " error = ", grad1 - grad2
END SUBROUTINE

END PROGRAM

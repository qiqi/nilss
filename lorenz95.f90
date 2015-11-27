PROGRAM Lorenz95_NILSS

    IMPLICIT NONE

    REAL(8), PARAMETER :: PI  = 4 * atan(1.0)
    REAL(8), PARAMETER :: s(1) = (/28.0/)
    INTEGER, PARAMETER :: nState = 40
    INTEGER, PARAMETER :: nChunks = 200, nStepsPerChunk = 500
    INTEGER, PARAMETER :: nHomoAdjoint = 20

    INTEGER :: iStep, iChunk, iAdjoint
    REAL(8) :: x(nState), y(nState, nHomoAdjoint+1), dotProductWeights(nState)
    REAL(8) :: grad(1, nHomoAdjoint+1), lss_grad(1)
    REAL(8), ALLOCATABLE :: history(:,:,:)

    INTEGER :: nTotalSteps
    REAL(8) :: windowFunction, timeFraction

    WRITE (*,*) 'TESTING ADJOINT IMPLEMENTATION'
    DO iChunk = 1, 10
        CALL testAdjoint
    END DO

    WRITE (*,*) 'INITIAL TRANSIENT'
    x(:) = 0.0
    x(1) = 1.0
    DO iStep = 1, 1000
        CALL Step(x, s)
    END DO

    WRITE (*,*) 'PRIMAL CALCULATION'
    Allocate(history(nState, nStepsPerChunk, nChunks))
    DO iChunk = 1, nChunks
        DO iStep = 1, nStepsPerChunk
            history(:, iStep, iChunk) = x(:)
            CALL Step(x, s)
        END DO
    END DO

    WRITE (*,*) 'ADJOINT CALCULATION'
    dotProductWeights(:) = 1.0
    ! NiLSS_init(nHomoAdjoint, nStateVariables,
    !            nDesignVariables, dotProductWeights)
    CALL NiLSS_init(nHomoAdjoint, nState, 1, dotProductWeights)

    y(:,:) = 0.0
    CALL RANDOM_NUMBER(y(:,1:nHomoAdjoint))

    DO iChunk = nChunks, 1, -1
        grad(:,:) = 0.0
        DO iStep = nStepsPerChunk, 1, -1
            x(:) = history(:, iStep, iChunk)
            DO iAdjoint = 1, nHomoAdjoint + 1
                CALL gradContribution(x, y(:,iAdjoint), s, grad(:,iAdjoint))
                CALL Adjoint(x, y(:,iAdjoint), s)
            END DO
            ! windowing
            nTotalSteps = nChunks * nStepsPerChunk
            timeFraction = real(iStep - 1 + (iChunk - 1) * nStepsPerChunk) &
                         / nTotalSteps
            windowFunction = sin(timeFraction * PI) * sin(timeFraction * PI) * 2
            y(nState, nHomoAdjoint + 1) = y(nState, nHomoAdjoint + 1) + &
                windowFunction / nTotalSteps
        END DO
        ! checkpoint the adjoint solution
        CALL NiLSS_checkpoint(y, grad)
    END DO

    CALL NiLSS_gradient(lss_grad)
    Write (*,*) 'NI-LSS Gradient = ', lss_grad(1)

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

SUBROUTINE Adjoint(x, y, s)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x(nState)
    REAL(8), INTENT(inout) :: y(nState)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.001
    INTEGER :: i, ip, ipp, im, imm

    REAL(8) :: dy(nState)
    DO i = 1, nState
        ip = mod(i, nState) + 1
        ipp = mod(i + 1, nState) + 1
        im = mod(i + nState - 2, nState) + 1
        imm = mod(i + nState - 3, nState) + 1
        dy(i) = -x(im) * y(ip) - x(ip) * y(ipp) &
              + x(imm) * y(im) + x(ipp) * y(ip) - y(i)
    END DO
    y(:) = y(:) + dt * dy(:)
END SUBROUTINE

SUBROUTINE gradContribution(x, y, s, grad)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x(nState), y(nState)
    REAL(8), INTENT(in) :: s(1)
    REAL(8), INTENT(inout) :: grad(1)

    REAL(8), PARAMETER :: dt = 0.001

    grad(1) = grad(1) + sum(y) * dt
END SUBROUTINE

SUBROUTINE testAdjoint
    IMPLICIT NONE
    REAL(8) :: x(nState,100), xp(nState), xm(nState), y(nState)
    REAL(8) :: s(1), sp(1), sm(1), eps, grad1, grad2(1)
    INTEGER :: iStep

    CALL RANDOM_NUMBER(x(:,1))
    CALL RANDOM_NUMBER(y)
    CALL RANDOM_NUMBER(s)
    x(:,1) = x(:,1) * 50
    s(:) = s(:) * 100

    eps = 1E-8

    xp(:) = x(:,1) + eps / 2
    xm(:) = x(:,1) - eps / 2

    CALL Step(xp, s)
    CALL Step(xm, s)
    grad1 = SUM((xp - xm) * y) / eps

    CALL Adjoint(x(:,1), y, s)
    grad2 = SUM(y)

    write (*,*) "FD=", grad1, "ADJ=", grad2, "ERR=", grad1 - grad2

    sp(:) = s(:) + eps / 2
    sm(:) = s(:) - eps / 2

    xp(:) = x(:,1)
    xm(:) = x(:,1)
    DO iStep = 2,100
        x(:,iStep) = x(:,iStep-1)
        CALL Step(x(:,iStep), s)
        CALL Step(xp, sp)
        CALL Step(xm, sm)
    END DO
    grad1 = SUM((xp - xm) * y) / eps

    grad2(1) = 0
    DO iStep = 99,1,-1
        CALL gradContribution(x(:,iStep), y, s, grad2)
        CALL Adjoint(x(:,iStep), y, s)
    END DO

    write (*,*) "FD=", grad1, "ADJ=", grad2, "ERR=", grad1 - grad2
END SUBROUTINE

END PROGRAM

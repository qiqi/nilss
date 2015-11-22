MODULE NiLSS

    REAL(8), ALLOCATABLE :: NiLSS_grad(:)

    REAL(8), ALLOCATABLE :: NiLSS_dotProductWeights(:)

    CONTAINS

    SUBROUTINE NiLSS_init(nHomoAdjoint, nStateVariables, &
                          nDesignVariables, dotProductWeights)
        INTEGER, INTENT(IN) :: nHomoAdjoint, nStateVariables, nDesignVariables
        REAL(8), INTENT(IN) :: dotProductWeights(:)

        ALLOCATE(NiLSS_grad(nDesignVariables))
    END SUBROUTINE

    SUBROUTINE NiLSS_Checkpoint(y, grad)
        REAL(8), INTENT(IN) :: y(:,:), grad(:,:)
    END SUBROUTINE

END MODULE NiLSS

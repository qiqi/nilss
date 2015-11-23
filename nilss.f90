MODULE NiLSS

    CONTAINS

    SUBROUTINE NiLSS_init(nHomoAdjoint, nStateVariables, &
                          nDesignVariables, dotProductWeights)
        INTEGER, INTENT(IN) :: nHomoAdjoint, nStateVariables, nDesignVariables
        REAL(8), INTENT(IN) :: dotProductWeights(:)
    END SUBROUTINE

    SUBROUTINE NiLSS_checkpoint(y, grad)
        REAL(8), INTENT(IN) :: y(:,:), grad(:,:)
    END SUBROUTINE

    SUBROUTINE NiLSS_gradient(grad)
        REAL(8), INTENT(OUT) :: grad(:)
    END SUBROUTINE

END MODULE NiLSS

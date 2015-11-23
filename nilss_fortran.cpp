#include<cassert>
#include "nilss.h"

nilss::NILSS * p_nilss;

extern "C" {
void nilss_init_(int* nHomoAdjoint, int* nStateVaraibles,
        int* nDesignVariables, double* dotProductWeights)
{
    assert(p_nilss == nullptr);
    p_nilss = new nilss::NILSS(*nHomoAdjoint, *nStateVaraibles,
            *nDesignVariables, dotProductWeights);
}

void nilss_checkpoint_(double ** y, double ** grad)
{
    assert(p_nilss != nullptr);
    p_nilss->checkpoint(y, grad);
}

void nilss_gradient_(double * grad)
{
    assert(p_nilss != nullptr);
    p_nilss->gradient(grad);
}
}

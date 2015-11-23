#include "nilss.h"

using namespace nilss;

NILSS::NILSS(int nHomoAdjoint, int nStateVariables,
                    int nDesignVariables, double * dotProductWeights)
{
}

NILSS::~NILSS()
{
}

void NILSS::checkpoint(double ** y, double ** grad)
{}

void NILSS::gradient(double * grad)
{}

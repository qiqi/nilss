#include <cstdlib>
#include <cxxtest/TestSuite.h>
#include "nilss.h"

using namespace nilss;

class TestNiLSS : public CxxTest::TestSuite
{
public:
    void testOrthogonalization(void)
    {
        double weights[] = {1.0, 1.0};
        NILSS lss(1, 2, 1, weights);

        double y[4];
        memset((void*)y, 0, 4 * sizeof(double));
        y[0] = 1.0;
        y[2] = 1.0;
        y[3] = 1.0;
        double grad[2];
        lss.checkpoint(y, grad);
        
        TS_ASSERT_DELTA(y[0], 1.0, 1E-12);
        TS_ASSERT_DELTA(y[1], 0.0, 1E-12);
        TS_ASSERT_DELTA(y[2], 0.0, 1E-12);
        TS_ASSERT_DELTA(y[3], 1.0, 1E-12);
    }
};

#include <cstdlib>
#include <iostream>
#include <cxxtest/TestSuite.h>
#include "nilss_solver.h"

using namespace Eigen;
using namespace nilss;

class TestNiLSS_Solver : public CxxTest::TestSuite
{
    private:

    void fillRDbc(int n, int m,
            std::vector<MatrixXd>& R,
            std::vector<MatrixXd>& D,
            std::vector<VectorXd>& b,
            std::vector<VectorXd>& c)
    {
        for (int i = 0; i < n; ++ i) {
            R.emplace_back(m,m);
            R.back().setIdentity();
            b.emplace_back(m);
            b.back().setOnes();
        }

        for (int i = 0; i <= n; ++ i) {
            D.emplace_back(m,m);
            D.back().setIdentity();
            c.emplace_back(m);
            c.back().setZero();
        }
    }

    public:

    void testNeutral(void)
    {
        int n = 100;
        std::vector<MatrixXd> R, D;
        std::vector<VectorXd> a, b, c;

        fillRDbc(n, 1, R, D, b, c);
        nilss_solve(R, D, b, c, a);

        TS_ASSERT_EQUALS(a.size(), n + 1);
        for (int i = 0; i <= n; ++ i) {
            TS_ASSERT_DELTA(a[i](0), i - n * 0.5, 1E-10);
        }
    }

    void testDecaying(void)
    {
        int n = 100;
        std::vector<MatrixXd> R, D;
        std::vector<VectorXd> a, b, c;

        fillRDbc(n, 1, R, D, b, c);

        double r = 0.01;
        for (int i = 0; i < n; ++ i) {
            R[i] *= r;
        }

        nilss_solve(R, D, b, c, a);

        TS_ASSERT_EQUALS(a.size(), n + 1);
        for (int i = 10; i <= n; ++ i) {
            TS_ASSERT_DELTA(a[i](0), 1.0 / (1 - r), 1E-10);
        }
    }

    void testGrowing(void)
    {
        int n = 100;
        std::vector<MatrixXd> R, D;
        std::vector<VectorXd> a, b, c;

        fillRDbc(n, 1, R, D, b, c);

        double r = 100.0;
        for (int i = 0; i < n; ++ i) {
            R[i] *= r;
        }

        nilss_solve(R, D, b, c, a);

        TS_ASSERT_EQUALS(a.size(), n + 1);
        for (int i = 0; i < n - 10; ++ i) {
            TS_ASSERT_DELTA(a[i](0), 1.0 / (1 - r), 1E-10);
        }
    }
};

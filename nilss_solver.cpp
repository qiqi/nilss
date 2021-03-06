#include <memory>
#include <cassert>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "nilss_solver.h"

using namespace Eigen;

namespace nilss {

std::unique_ptr<MatrixXd> assemble_kkt(
        const std::vector<MatrixXd>& R, const std::vector<MatrixXd>& D)
{
    int n = R.size(), m = R[0].rows();
    assert(n > 0 && D.size() == n + 1);

    int kktSize = (2 * n + 1) * m;
    MatrixXd * pMat = new MatrixXd(kktSize, kktSize);
    pMat->setZero();

    for (int i = 0; i <= n; ++ i) {
        pMat->block(i * m, i * m, m, m) = D[i];
    }

    int halfSize = (n + 1) * m;
    for (int i = 0; i < n; ++ i) {
        pMat->block(halfSize + i * m, (i + 1) * m, m, m) = MatrixXd::Identity(m, m);
        pMat->block(halfSize + i * m, i * m, m, m) = -R[i];

        pMat->block((i + 1) * m, halfSize + i * m, m, m) = MatrixXd::Identity(m, m);
        pMat->block(i * m, halfSize + i * m, m, m) = -R[i].transpose();
    }

    return std::unique_ptr<MatrixXd>(pMat);
}

std::unique_ptr<VectorXd> assemble_rhs(
        const std::vector<VectorXd>& b, const std::vector<VectorXd>& c)
{
    int n = b.size(), m = b[0].size();
    assert(n > 0 && c.size() == n + 1);

    int kktSize = (2 * n + 1) * m;
    VectorXd * pVec = new VectorXd(kktSize);

    for (int i = 0; i <= n; ++ i) {
        pVec->segment(i * m, m) = c[i];
    }

    int halfSize = (n + 1) * m;
    for (int i = 0; i < n; ++ i) {
        pVec->segment(halfSize + i * m, m) = b[i];
    }

    return std::unique_ptr<VectorXd>(pVec);
}

void nilss_solve(const std::vector<MatrixXd>& R,
                 const std::vector<MatrixXd>& D,
                 const std::vector<VectorXd>& b,
                 const std::vector<VectorXd>& c,
                 std::vector<VectorXd>& a)
{
    int n = R.size();
    assert(D.size() == n + 1);
    assert(b.size() == n);
    assert(c.size() == n + 1);

    std::unique_ptr<MatrixXd> kkt = assemble_kkt(R, D);
    std::unique_ptr<VectorXd> rhs = assemble_rhs(b, c);

    typedef SparseMatrix<double> SpMat;
    SpMat A(kkt->sparseView());

    SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(A); 
    solver.factorize(A); 
    VectorXd sol = solver.solve(*rhs); 

    //VectorXd sol = kkt->partialPivLu().solve(*rhs);

    assert(sol.size() % (2 * n + 1) == 0);
    int m = sol.size() / (2 * n + 1);

    a.empty();
    a.reserve(n + 1);
    for (int i = 0; i <= n; ++ i) {
        a.push_back(sol.segment(i * m, m));
    }
}

}

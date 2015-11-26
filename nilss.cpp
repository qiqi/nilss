#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "nilss.h"
#include "nilss_solver.h"

using namespace nilss;

NILSS::NILSS(int nHomoAdjoint, int nStateVariables,
             int nDesignVariables, const double * dotProductWeights)
    : nHomo_(nHomoAdjoint), size_(nStateVariables), nGrad_(nDesignVariables)
{
    for (int i = 0; i < size_; ++i) {
        dotWeights_.push_back(dotProductWeights[i]);
    }
}

double NILSS::dotProd_(const double * y1, const double * y2) const
{
    double dot = 0;
    for (int i = 0; i < size_; ++ i) {
        dot += dotWeights_[i] * y1[i] * y2[i];
    }
    return dot;
}

void NILSS::axpy_(double * y, const double * x, double a) const
{
    for (int i = 0; i < size_; ++ i) {
        y[i] += a * x[i];
    }
}

void NILSS::scale_(double * y, double a) const
{
    for (int i = 0; i < size_; ++ i) {
        y[i] *= a;
    }
}

void NILSS::checkpoint(double * y, const double * grad)
{
    // Gram Schmidt orthonormalization
    R_.emplace_back(nHomo_, nHomo_);
    Eigen::MatrixXd & R = R_.back();
    for (int i = 0; i < nHomo_; ++ i) {
        R.col(i).setZero();
        double * yi = y + i * size_;
        for (int j = 0; j < i; ++j) {
            double * yj = y + j * size_;
            R(j,i) = dotProd_(yi, yj);
            axpy_(yi, yj, -R(j,i));
        }
        R(i,i) = sqrt(dotProd_(yi, yi));
        scale_(yi, 1.0/R(i,i));
    }
    // Orthogonalization
    b_.emplace_back(nHomo_);
    Eigen::VectorXd & b = b_.back();
    b.setZero();
    double * yInHomo = y + nHomo_ * size_;
    for (int j = 0; j < nHomo_; ++j) {
        double * yj = y + j * size_;
        b(j) = dotProd_(yj, yInHomo);
        axpy_(yInHomo, yj, -b(j));
    }
    // Store gradient
    stored_grad_.emplace_back(nHomo_ + 1, nGrad_);  // homo and inhomo
    Eigen::MatrixXd & stored_grad = stored_grad_.back();
    for (int i = 0; i <= nHomo_; ++i) {
        for (int j = 0; j < nGrad_; ++j) {
            stored_grad(i,j) = grad[i * nGrad_ + j];
        }
    }
}

double NILSS::window_(double x) const
{
    const double PI = atan(1.0) * 4;
    double s = sin(x * PI);
    return s * s;
}

void NILSS::gradient(double * gradient) const
{
    std::vector<Eigen::MatrixXd> identities;
    std::vector<Eigen::VectorXd> zeros;
    identities.reserve(R_.size() + 1);
    zeros.reserve(R_.size() + 1);

    for (size_t i = 0; i <= R_.size(); ++ i) {
        double window = window_((double)i / R_.size());
        identities.emplace_back(nHomo_, nHomo_);
        identities.back().setIdentity();
        identities.back() *= window;
        zeros.emplace_back(nHomo_);
    }

    std::vector<Eigen::VectorXd> a;
    nilss_solve(R_, identities, b_, zeros, a);
    assert(a.size() == stored_grad_.size());

    Eigen::VectorXd window(a.size());
    for (size_t i = 0; i < a.size(); ++ i) {
        window(i) = window_((double)i / (R_.size() - 2));
    }
    window /= window.sum();

    // combine the gradients
    Eigen::VectorXd grad(nGrad_);
    grad.setZero();
    for (size_t i = 0; i < a.size(); ++ i) {
        auto gradi = stored_grad_[i];
        for (int j = 0; j <= nHomo_; ++ j) {
            double aij = (j < nHomo_) ? a[i](j) : 1.0;
            grad += aij * gradi.row(j);
        }
    }

    for (int i = 0; i < nGrad_; ++i) {
        gradient[i] = grad(i);
    }
}

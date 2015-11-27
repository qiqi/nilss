#ifndef NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_SOLVER_H
#define NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_SOLVER_H

#include<Eigen/Dense>

namespace nilss
{
    void nilss_solve(const std::vector<Eigen::MatrixXd>& R,
                     const std::vector<Eigen::MatrixXd>& D,
                     const std::vector<Eigen::VectorXd>& b,
                     const std::vector<Eigen::VectorXd>& c,
                     std::vector<Eigen::VectorXd>& a);
}

#endif

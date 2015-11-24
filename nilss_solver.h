#ifndef NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_SOLVER_H
#define NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_SOLVER_H

#include<Eigen/Dense>

void nilss_solve(const std::vector<Eigen::MatrixXd>& R_,
                 const std::vector<Eigen::MatrixXd>& identities,
                 const std::vector<Eigen::VectorXd>& b_,
                 const std::vector<Eigen::VectorXd>& zeros,
                 std::vector<Eigen::VectorXd>& a);

#endif

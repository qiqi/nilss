#ifndef NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_H
#define NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_H

#include <Eigen/Dense>

namespace nilss {
    class NILSS {
        private:
            int nHomo_;
            int size_;
            int nGrad_;
            std::vector<double> dotWeights_;
            std::vector<Eigen::MatrixXd> R_;
            std::vector<Eigen::VectorXd> b_;
            std::vector<Eigen::MatrixXd> stored_grad_;

            double dotProd_(const double * y1, const double * y2) const;
            void axpy_(double * y, const double * x, double a) const;
            double window_(double x) const;

        public:
            NILSS(int nHomoAdjoint, int nStateVariables,
                  int nDesignVariables, const double * dotProductWeights);

            void checkpoint(double* const* y, const double* const* grad);

            void gradient(double * grad) const;
    };
}

#endif

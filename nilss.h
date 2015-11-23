#ifndef NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_H
#define NON_INTRUSIVE_LEAST_SQUARES_SHADOWING_NILSS_H

namespace nilss {
    class NILSS {
        public:
            NILSS(int nHomoAdjoint, int nStateVariables,
                    int nDesignVariables, double * dotProductWeights);

            virtual ~NILSS();

            void checkpoint(double ** y, double ** grad);

            void gradient(double * grad);
    };
}

#endif

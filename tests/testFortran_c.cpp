#include <iostream>
#include <cassert>

extern "C" {
void c_subroutine_(double * x, double * y, double * z)
{
    std::cout << "Testing 1D array x" << std::endl;
    assert(x[0] == 1.0);
    assert(x[1] == 2.0);

    for (int i = 0; i < 3; ++i) {
        std::cout << "Testing 2D array y(" << i << ",:)" << std::endl;
        assert(y[i]   ==  1 + i);
        assert(y[3+i] == 11 + i);

        std::cout << "Testing 2D allocatable array z(" << i << ",:)" << std::endl;
        assert(z[i]   ==  1 + i);
        assert(z[3+i] == 11 + i);
    }
}
}

#include <iostream>
#include "CSR_matrix.h"
#include <cmath>


int main()
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = Simple_Iteration_Method(A, b, x0, tol, 10000);
    std::cout << x;
}


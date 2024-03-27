#include "../my/CSR_matrix.h"

#include <gtest/gtest.h>


TEST(CSR_matrix, constructor_matrix)
{
    vec<int> valuess = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(valuess , 3);
    CSR_matrix<int> CSR(M);

    vec<int> vals = {1, 3, 5, 7, 9};
    vec<std::size_t> cols = {0, 2, 1, 0, 2};
    vec<std::size_t> rows = {0, 2, 3, 5};

    EXPECT_EQ(true, CSR.check(vals, cols, rows));
}

TEST(CSR_matrix, get_element)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_EQ(CSR(i, j), vals[i * 3 + j]);
}

TEST(CSR_matrix, get_cols)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    vec<int> expected = {0, 2, 1, 0, 2}; 
    EXPECT_EQ(CSR(0, 0), 1);
    EXPECT_EQ(CSR(0, 1), 0);
}

TEST(CSR_matrix, mul)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);

    vec<int> x = {1, 2, 3};
    x = CSR * x;
    vec<int> expected = {14, 32, 50};
    for (std::size_t j = 0; j < 3; ++j)
            EXPECT_EQ(x[j], expected[j]);
}

TEST(CSR_matrix, Simple_Iteration_Method)
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = Simple_Iteration_Method(A, b, x0, tol, 10000);
    vec<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR_matrix, SIM_Chebyshev)
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = SIM_Chebyshev(A, b, x0, tol, 10000, 3, 0.287, 10.261);
    vec<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR_matrix, Jacobi_Method)
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = Jacobi_Method(A, b, x0, tol, 10000);
    vec<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR_matrix, Gauss_Zeidel_Method)
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = Gauss_Zejdel_Method(A, b, x0, tol, 10000);
    vec<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR_matrix , lambda_max)
{
    vec<double> vals = {76, 0, 2, 0, 3, 58, 0, 0, 0, 0, 84, 0, 3, 0, 0, 64}; 
    Matrix<double> M(vals , 4);
    CSR_matrix<double> A(M);
    EXPECT_NEAR(A.lambda_max(), 84, 1e-6);
}

TEST(CSR_matrix, Steepest_Descent)
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = Steepest_Descent(A, b, x0, tol, 10000);
    vec<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], x[j], 0.01);
}

TEST(CSR_matrix, Sym_Gauss_Zeidel_Method)
{
    vec<double> vals = {1, 2, 0, 2, 6, 1, 0, 1, 10};
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {4, 4, 4};
    vec<double> b = {5, 17, 32};
    double tol = 1e-20;
    vec<double> x = Sym_Gauss_Zejdel_Method(A, b, x0, tol, 10000);
    vec<double> expected = {1, 2, 3};
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], x[j], 0.01);
}
/*
TEST(CSR_matrix, kr5)
{
    vec<double> vals = {50, 0, 0, 6, 0, 87, 6, 0, 0, 0, 59, 0, 0, 0, 6, 74}; 
    Matrix<double> M(vals , 3);
    CSR_matrix<double> A(M);

    vec<double> x0 = {0, 0, 0, 0};
    vec<double> b = {1, 1, 1, 1};
    double tol = 1e-12;

    vec<std::size_t> it = Simple_Iteration_Method(A, b, x0, tol, 10000);

    
}



TEST(CSR_matrix, kr1)
{
    vec<double> vals = {0, 0, 0, 2, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
    Matrix<double> M(vals , 4);
    CSR_matrix<double> A(M);

    std::cout << A.get_vals() << std::endl;
    std::cout << A.get_cols() << std::endl;
    std::cout << A.get_rows();
}


TEST(CSR_matrix, kr3)
{
    vec<double> vals = {79, 9, 0, 6, 0, 55, 2, 6, 0, 3, 72, 3, 8, 0, 0, 81};
    Matrix<double> M(vals , 4);
    CSR_matrix<double> A(M);

    vec<double> x0 = {367326700, 0, 0, 0};
    vec<double> b = {5, 3, 2, 7};
    double tol = 1e-20;
    std::size_t it = Gauss_Zejdel_Method_kr(A, b, x0, tol, 10000);
    std::cout << it << std::endl;
}

TEST(CSR_matrix, get_rows)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    vec<int> expected = {0, 2, 3, 5};
    EXPECT_EQ(true, CSR.get_rows() == expected[0]);
}*/

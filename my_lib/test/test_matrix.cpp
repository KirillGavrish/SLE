//#include "../my/vec.h"
#include "../my/Matrix.h"

#include <gtest/gtest.h>


TEST(Matrix, operator_get_element)
{
    Matrix<int> test_matrix({1, 2, 3, 4}, 2);
    EXPECT_EQ(test_matrix(0, 0), 1);
    EXPECT_EQ(test_matrix(1, 1), 4);
    EXPECT_EQ(test_matrix(0, 1), 2);
    EXPECT_EQ(test_matrix(1, 0), 3);
}

TEST(Matrix, get_vals)
{
    vec<int> vals = {1, 2, 3, 4};
    Matrix<int> test_matrix(vals, 2);
    EXPECT_EQ(test_matrix.vals(), vals);
}
/*
TEST(Matrix, get_pos_vals)
{
    vec<int> vals = {-1, 2, -3, 4};
    Matrix<int> test_matrix(vals, 2);
    vec<int> pos_vals = {2, 4};
    EXPECT_EQ(test_matrix.pos_vals(), pos_vals);
}
*/
TEST(Matrix, constructor_m_n_default)
{
    Matrix<int> M;
    EXPECT_EQ(M.vals().size(), 0);
}

TEST(Matrix, constructor_m_n_zero)
{
    Matrix<int> M(0 , 0);
    EXPECT_EQ(M.vals().size(), 0);
    EXPECT_EQ(M.get_width(), 0);
}

TEST(Matrix, constructor_m_n)
{
    Matrix<int> M(2, 2);
    EXPECT_EQ(M.vals().size(), 4);
    EXPECT_EQ(M.get_width(), 2);
}

TEST(Matrix, constructor_vec)
{
    vec<int> vals = {-1, 2, -3, 4};
    Matrix<int> M(vals , 2);
    EXPECT_EQ(M(0, 1), 2);
}

TEST(Matrix, get_width)
{
    vec<int> vals = {-1, 2, -3, 4};
    Matrix<int> M(vals , 2);
    EXPECT_EQ(M.get_width(), 2);
}

TEST(Matrix, get_height)
{
    vec<int> vals = {-1, 2, -3, 4};
    Matrix<int> M(vals , 2);
    EXPECT_EQ(M.get_height(), 2);
}

TEST(Matrix, mul_vec)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M(vals , 3);
    vec<int> v = {1, 2, 3};
    vec<int> product = M * v;
    vec<int> expected = {14, 32, 50};
    EXPECT_EQ(product, expected);
}

TEST(Matrix, mul_Matrix)
{
    Matrix<int> mtr1 = Matrix<int>({4, 5, 6, 7, 8, 10, 12, 14, 12, 15, 18, 21}, 4);
	Matrix<int> mtr2 = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, 3);
	Matrix<int> res = mtr1 * mtr2;
    Matrix<int> mtr = Matrix<int>({136, 158, 180, 272, 316, 360, 408, 474, 540}, 3);
    for (std::size_t i = 0; i < res.get_width(); ++i)
            EXPECT_EQ(mtr.get_col(i), res.get_col(i));
}

TEST(Matrix, transpose)
{
	Matrix<int> mtr = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, 3);
	Matrix<int> res = Matrix<int>({1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12}, 4);
    mtr = mtr.transpose();
	for (std::size_t i = 0; i < res.get_width(); ++i)
            EXPECT_EQ(mtr.get_col(i), res.get_col(i));
}
 
TEST(Matrix, QR_decomp)
{
    Matrix<double> mtr = Matrix<double>({1, 3, 5, 1, 3, 1, 2, -1, 7}, 3);
    std::pair<Matrix<double>, Matrix<double>> QR = QR_decomp(mtr);
    Matrix<double> res = QR.first * QR.second;
    vec<double> v1 = {-0.408248290463863, -0.577350269189626, 0.707106781186548, -0.408248290463863, -0.577350269189626, -0.707106781186548, -0.816496580927726, 0.577350269189626, 0};
    vec<double> v2 = {-2.449489742783178, -1.632993161855452, -8.16496580927726, 0, -4.04145188432738, 0.577350269189626, 0, 0, 2.82842712474619};
    for(size_t i = 0; i < v1.size(); i++){
        ASSERT_NEAR(QR.first.vals()[i], v1[i], std::abs(0.001 * QR.first.vals()[0]));
        ASSERT_NEAR(QR.second.vals()[i], v2[i], std::abs(0.001 * QR.second.vals()[0]));
        ASSERT_NEAR(res.vals()[i], mtr.vals()[i], std::abs(0.001 * mtr.vals()[0]));
    }
}

TEST(Matrix, add_col_H)
{
    Matrix<int> mtr = Matrix<int>({1, 2, 3, 4}, 2);
    Matrix<int> res = Matrix<int>({1, 2, 1, 3, 4, 2, 0, 0, 3}, 3);
    mtr = add_col_H(mtr, {1, 2, 3});
    for (std::size_t i = 0; i < res.get_width(); ++i)
            EXPECT_EQ(mtr.get_col(i), res.get_col(i));
}

TEST(Matrix, add_col)
{
    Matrix<int> mtr = Matrix<int>({1, 2, 3, 4}, 2);
    Matrix<int> res = Matrix<int>({1, 2, 1, 3, 4, 2}, 3);
    mtr = add_col(mtr, {1, 2});
    for (std::size_t i = 0; i < res.get_width(); ++i)
            EXPECT_EQ(mtr.get_col(i), res.get_col(i));
}

TEST(Matrix, get_col)
{
    Matrix<int> mtr = Matrix<int>({1, 2, 3, 4}, 2);
    Matrix<int> res = Matrix<int>({1, 2, 1, 3, 4, 2}, 3);
    mtr = add_col(mtr, res.get_col(2));
    for (std::size_t i = 0; i < res.get_width(); ++i)
            EXPECT_EQ(mtr.get_col(i), res.get_col(i));
}

TEST(Matrix, Inverse_Gauss_Method)
{
    Matrix<double> M = Matrix<double>({1, 2, 3, 0, 5, 6, 0, 0, 9}, 3);
    vec<double> b = {35, 67, 63};
    vec<double> expected = {4, 5, 7};
    vec<double> res = Inverse_Gauss_Method(M, b);
    for (std::size_t j = 0; j < 3; ++j)
        EXPECT_NEAR(expected[j], res[j], 0.001);
}
/*
TEST(Matrix, mul_other)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M(vals , 3);
    M *= M;
    vec<int> expected_vals = {30, 36, 42, 66, 81, 96, 102, 126, 150};
    EXPECT_EQ(M.vals(), expected_vals);
    EXPECT_EQ(M.get_width(), 3);
}

TEST(Matrix, add_other)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M(vals , 3);
    M = M + M;
    vec<int> expected_vals = {2, 4, 6, 8, 10, 12, 14, 16, 18};
    EXPECT_EQ(M.vals(), expected_vals);
}

TEST(Matrix, sub_other)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M1(vals , 3);
    Matrix<int> M2(vals , 3);
    M1 -= M2;
    vec<int> expected_vals = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    EXPECT_EQ(M1.vals(), expected_vals);
}

TEST(Matrix, mul_other_2)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M(vals , 3);
    M = M * M;
    vec<int> expected_vals = {30, 36, 42, 66, 81, 96, 102, 126, 150};
    EXPECT_EQ(M.vals(), expected_vals);
    EXPECT_EQ(M.get_width(), 3);
}

TEST(Matrix, add_other_2)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M(vals , 3);
    M = M + M;
    vec<int> expected_vals = {2, 4, 6, 8, 10, 12, 14, 16, 18};
    EXPECT_EQ(M.vals(), expected_vals);
}

TEST(Matrix, sub_other_2)
{
    vec<int> vals = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix<int> M1(vals , 3);
    Matrix<int> M2(vals , 3);
    M1 = M1 - M2;
    vec<int> expected_vals = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    EXPECT_EQ(M1.vals(), expected_vals);
} */

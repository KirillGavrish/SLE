//#include "../my/vec.h"
#include "../my/Matrix.h"

#include <gtest/gtest.h>


TEST(Matrix, operator_get_element)
{
    Matrix<int> test_matrix({1, 2, 3, 4}, 2);
    EXPECT_EQ(test_matrix(0, 0), 1);
    EXPECT_EQ(test_matrix(1, 1), 4);
}

TEST(Matrix, get_vals)
{
    vec<int> vals = {1, 2, 3, 4};
    Matrix<int> test_matrix(vals, 2);
    EXPECT_EQ(test_matrix.vals(), vals);
}

TEST(Matrix, get_pos_vals)
{
    vec<int> vals = {-1, 2, -3, 4};
    Matrix<int> test_matrix(vals, 2);
    vec<int> pos_vals = {2, 4};
    EXPECT_EQ(test_matrix.pos_vals(), pos_vals);
}

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
    vec<int> product = M * vals;
    vec<int> expected = {14, 32, 50};
    EXPECT_EQ(product, expected);
}

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
    M += M;
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
}

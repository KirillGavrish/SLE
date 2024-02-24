#include "../my/CSR_matrix.h"

#include <gtest/gtest.h>

/*
TEST(CSR_matrix, constructor_matrix)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    CSR_matrix<int> expected;
    expected = {{1, 3, 5, 7, 9}, {0, 2, 1, 0, 2}, {0, 2, 3, 5}};
    EXPECT_EQ(CSR, expected);
}
*/

TEST(CSR_matrix, get_vals)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    EXPECT_EQ(CSR.get_vals(), M.pos_vals());
}

TEST(CSR_matrix, get_element)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    EXPECT_EQ(CSR(0, 0), 1);
    EXPECT_EQ(CSR(0, 1), 0);
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

TEST(CSR_matrix, get_rows)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    vec<int> expected = {0, 2, 3, 5};
    EXPECT_EQ(CSR(0, 0), 1);
    EXPECT_EQ(CSR(0, 1), 0);
}

TEST(CSR_matrix, insert_check)
{
    vec<int> vals = {1, 0, 3, 0, 5, 0, 7, 0, 9};
    Matrix<int> M(vals , 3);
    CSR_matrix<int> CSR(M);
    std::size_t b = 2;
    CSR.insert(0, 1, b);
    EXPECT_EQ(CSR(0, 1), b);
    EXPECT_EQ(CSR.get_vals().size(), 6);
    EXPECT_EQ(CSR.get_rows()[1], 3);
    EXPECT_EQ(CSR.get_rows()[2], 4);
    EXPECT_EQ(CSR.get_rows()[1], 6);
}

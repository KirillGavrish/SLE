#include <gtest/gtest.h>

#include "../my/vec.h"
//#include <my/Matrix.h>
//#include <my/CSR_matrix.h>


TEST(vec, add) 
{
    vec<int> a = {1, 4};
    vec<int> b = {5, 2};
    vec<int> c = a + b;
    EXPECT_EQ(c[0], 6);
    EXPECT_EQ(c[1], 6);
}

TEST(vec, mul_right)
{
    vec<int> a = {4, 1};
    vec<int> b =  a * 5;
    vec<int> c = {20, 5};
    EXPECT_EQ(b[0], 20);
    EXPECT_EQ(b[1], 5);
}

TEST(vec, sub)       
{
    vec<int> a = {1, 4};
    vec<int> b = {5, 2};
    vec<int> c = a + b;
    EXPECT_EQ(c[0], 6);
    EXPECT_EQ(c[1], 6);
}

TEST(vec, mul_left)
{
    vec<int> a = {4, 1};
    vec<int> b = 5 * a;
    vec<int> c = {20, 5};
    EXPECT_EQ(b[0], 20);
    EXPECT_EQ(b[1], 5);
}

TEST(vec, dot)
{
    vec<int> a = {1, 4};
    vec<int> b = {5, 2};
    int c = dot(a, b);
    EXPECT_EQ(c, 13);
}

TEST(vec, max)
{
    vec<int> a = {1, 4};
    int c = max(a);
    EXPECT_EQ(c, 4);
}

#include <gtest/gtest.h>

#include "../my/vec.h"
//#include <my/Matrix.h>
//#include <my/CSR_matrix.h>


TEST(vec, plus) 
{
    vec<int> a = {1, 4};
    vec<int> b = {5, 2};
    vec<int> c = a + b;
    EXPECT_EQ(c[0], 6);
    EXPECT_EQ(c[1], 6);
}

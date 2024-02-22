#include <gtest/gtest.h>

#include <gtest/gtest.h>

#include <my/code/vec.h>
//#include <my/code/Matrix.h>
//#include <my/code/CSR_matrix.h>


TEST(vec, plus) 
{
    vec<int> a = {1};
    vec<int> b = {5};
    vec<int> c = a + b;
    EXPECT_IQ(c[0], 6);
}

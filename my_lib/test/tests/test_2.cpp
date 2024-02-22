#include <gtest/gtest.h>

#include <my/code/vec.h>
#include <my/code/Matrix.h>
#include <my/code/CSR_matrix.h>

TEST(vec, dot) {
    vec<int> a = {5, 5};
    vec<int> b = a;
    EXPECT_IQ(a dot b, 25);
}



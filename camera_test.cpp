#include "gtest/gtest.h"
#include "camera.h"
#include <iostream>

namespace idl{
namespace geometry{

TEST(EmptyLayer, Construction)
{
    double a[3] = {0.1, 0.2, 0.3};
    BrownRadialDistortiond lens(a);
    double p[3] = {1, 2, 3};
    double q[3], j1[9], j2[9];
    lens.project(p, q, j1, j2);
}

TEST(BrownRadialDistortiond, jacobRay)
{
    double a[3] = {0.1, 0.2, 0.3};
    BrownRadialDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    double q0[3], jRay[9];
    
    lens.project(p0, q0, NULL, jRay);
    double dx = 1e-6;
    for (int i = 0; i < 3; ++i)
    {
        double p1[3], q1[3];
        memcpy(p1, p0, 3 * sizeof(double));
        p1[i] += dx;
        lens.project(p1, q1);
        for (int j = 0; j < 3; ++j)
            ASSERT_NEAR((q1[j] - q0[j]) / dx, jRay[j + i * 3], 1e-5);
    }

}


} // geometry
} // idl

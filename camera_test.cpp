#include "gtest/gtest.h"
#include "camera.h"
#include <iostream>

namespace idl{
namespace geometry{

template<typename Distortion>
void testJacobRay(const Distortion& distortion, typename Traits<Distortion>::Scalar ray[3])
{
    typedef typename Traits<Distortion>::Scalar Scalar;

    Scalar refracted[3], jRay[9];
    distortion.project(ray, refracted, NULL, jRay);

    Scalar dx = 1e-6;
    for (int i = 0; i < 3; ++i)
    {
        Scalar p[3], q[3];
        memcpy(p, ray, 3 * sizeof(Scalar));
        p[i] += dx;
        distortion.project(p, q);
        for (int j = 0; j < 3; ++j)
        {
            ASSERT_NEAR((q[j] - refracted[j]) / dx, jRay[j + 3 * i], Scalar(1e-5));
        }
    }
}

template<typename Distortion>
void testJacobDistort(const Distortion& distortion, typename Traits<Distortion>::Scalar ray[3])
{
    typedef typename Traits<Distortion>::Scalar Scalar;
    enum { kDataSize = Traits<Distortion>::kDataSize };

    Scalar refracted[3], jDistort[3 * kDataSize];
    distortion.project(ray, refracted, jDistort, NULL);

    Scalar dx = 1e-6;
    std::cout << kDataSize << std::endl;
    for (int i = 0; i < kDataSize; ++i)
    {
        Scalar coeff[kDataSize], q[3];
        memcpy(coeff, distortion.data(), kDataSize * sizeof(Scalar));
        coeff[i] += dx;
        Distortion another(coeff);
        another.project(ray, q);
        for (int j = 0; j < 3; ++j)
            ASSERT_NEAR((q[j] - refracted[j]) / dx, jDistort[j + 3 * i], 1e-6);
            //std::cout << (q[j] - refracted[j]) / dx << " " << jDistort[j + i * 3] << std::endl;;
    }
}

TEST(EmptyLayer, Construction)
{
    double a[3] = {0.1, 0.2, 0.3};
    BrownDistortiond lens(a);
    double p[3] = {1, 2, 3};
    double q[3], j1[9], j2[9];
    lens.project(p, q, j1, j2);
}

TEST(BrownDistortiond, jacobRay)
{
    double a[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    BrownDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    testJacobRay(lens, p0);
}

TEST(BrownDistortiond, jacobDistort)
{
    double a[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    BrownDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    testJacobDistort(lens, p0);
}

TEST(CatadioptricDistortiond, project)
{
    double a[1] = {1};
    CatadioptricDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    double q0[3];
    lens.project(p0, q0);
    std::cout << q0[0] << " " << q0[1] << " " << q0[2] << std::endl;
}

TEST(CatadioptricDistortiond, jacobRay)
{
    double a[1] = {1};
    CatadioptricDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    testJacobRay(lens, p0);
}

TEST(CatadioptricDistortiond, jacobDistort)
{
    double a[1] = {2};
    CatadioptricDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    testJacobRay(lens, p0);
}

TEST(BrownCatadioptricDistortiond, jacobRay)
{
    double a[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    BrownCatadioptricDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    testJacobRay(lens, p0);
}

TEST(BrownCatadioptricDistortiond, jacobDistort)
{
    double a[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    BrownCatadioptricDistortiond lens(a);
    double p0[3] = {1, 2, 3};
    testJacobDistort(lens, p0);
}
} // geometry
} // idl

#include <gtest/gtest.h>
#include "pose.h"

namespace idl {
namespace geometry {

TEST(Pose, RotationMatrixConstruction)
{
    RotationMatrix<double, OwnedStorage> R0;

    double data[9] = {
        -0.6949,    0.7135,    0.0893,
        -0.1920,   -0.3038,    0.9332,
         0.6930,    0.6313,    0.3481};
    double data2[9] = {
        -6949,    7135,    893,
        -1920,   -3038,    9332,
         6930,    6313,    3481};
    RotationMatrix<double, OwnedStorage> R1(data);
    RotationMatrix<double, SharedStorage> R2(data2);

    R1 = R2;
    R2 = R1;
    RotationMatrix<double> R3(R2);
    R0 = R1;
}

TEST(Pose, RotationMatrixInversionCombination)
{
    double data[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };
    RotationMatrix<double> R1(data), R2 = R1.inverse();
    RotationMatrix<double> C = R1.combine(R2);
    ASSERT_NEAR(1, C.data()[0], 1e-3);
    ASSERT_NEAR(1, C.data()[4], 1e-3);
    ASSERT_NEAR(1, C.data()[8], 1e-3);
}

TEST(Pose, RotationMatrixTransformInOut)
{
    double p0[3] = {1, 2, 3}, p1[3], p2[3];
    double data[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };
    RotationMatrix<double> R(data);
    R.transformOut(p0, p1);
    R.transformOut(p1, p2);
    ASSERT_NEAR(p0[0], p1[0], 1e-3);
    ASSERT_NEAR(p0[1], p1[1], 1e-3);
    ASSERT_NEAR(p0[2], p1[2], 1e-3);
}

TEST(Pose, RotationMatrixTransformInJacob)
{
    double p0[3] = {1, 2, 3}, p1[3], p2[3];
    double jacob[3*9];
    double data[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };
    RotationMatrix<double> R0(data);
    R0.transformOut(p0, p1, jacob);
    double dx = 1e-5;
    for (int i = 0; i < 9; ++i)
    {
        RotationMatrix<double> R1(data);
        R1.data()[i] += dx;
        R1.transformOut(p0, p2);
        ASSERT_NEAR(jacob[i * 3 + 0], (p2[0] - p1[0]) / dx,1e-5);
        ASSERT_NEAR(jacob[i * 3 + 1], (p2[1] - p1[1]) / dx,1e-5);
        ASSERT_NEAR(jacob[i * 3 + 2], (p2[2] - p1[2]) / dx,1e-5);
    }
}

TEST(Pose, RotationMatrixRelate)
{
    double data[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };
    RotationMatrix<double> R0(data), R1(data);
    RotationMatrix<double> R0sqr = R0.combine(R0);
    RotationMatrix<double> R0sqr2 = R0.inverse().relate(R0);
    for (int i = 0; i < 9; ++i)
        ASSERT_DOUBLE_EQ(R0sqr.data()[i], R0sqr2.data()[i]);

}

TEST(Pose, AngleAxisConstruction)
{
    double data[3] = {1, 2, 3};
    AngleAxis<double> aa0(data);
    AngleAxis<double, SharedStorage> aa1(data);
    AngleAxis<double> aa2;
    aa1 = aa0;
}

TEST(Pose, AngleAxisToRotationMatrix)
{
    double groundtruth[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };

    double data[3] = {1, 2, 3};
    AngleAxis<double> aa(data);
    double rm[9];
    aa.toRotationMatrix(rm);
    for (int i = 0; i < 9; ++i)
        ASSERT_NEAR(groundtruth[i], rm[i], 1e-9);
}

TEST(Pose, AngleAxisToRotationMatrixJacob)
{
    double groundtruth[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };

    double data[3] = {1, 2, 3};
    AngleAxis<double> aa(data);
    double rm0[9], jacob0[27];
    aa.toRotationMatrix(rm0, jacob0);
    double dx = 1e-5;
    for (int i = 0; i < 3; ++i)
    {
        AngleAxis<double> aa(data);
        aa.data()[i] += dx;
        double rm1[9];
        aa.toRotationMatrix(rm1);
        for (int j = 0; j < 9; ++j)
            ASSERT_NEAR((rm1[j] - rm0[j]) / dx, jacob0[9 * i + j], 1e-5);
    }
}

TEST(Pose, RotationMatrixToAngleAxis)
{
    double data[9] = {
        -0.694920557641312,
         -0.192006972791999,
         0.692978167741770,
         0.713520990527788,
         -0.303785044339471,
         0.631349699383718,
         0.089292858861913,
         0.933192353823647,
         0.348107477830265 };
    RotationMatrix<double> R(data);
    double aa0[3], jacob0[27];
    R.toAngleAxis(aa0, jacob0);
    double dx = 1e-8;
    for (int i = 0; i < 9; ++i)
    {
        RotationMatrix<double> R(data);
        R.data()[i] += dx;
        double aa1[3];
        R.toAngleAxis(aa1);
        for (int j = 0; j < 3; ++j)
        {
            ASSERT_NEAR((aa1[j] - aa0[j]) / dx, jacob0[3 * i + j], 1e-5);
        }
    }
}

TEST(Pose, AngleAxisInversionCombination)
{
    AngleAxis<double> aa0(1, 2, 3);
    AngleAxis<double> aa1 = aa0.inverse();
    AngleAxis<double> C = aa0.combine(aa1);
}

TEST(Pose, PoseConstruction)
{
    AngleAxis<double> aa(1, 2, 3);
    double loc[3] = {4, 5, 6};
    Pose<double, AngleAxis> P(aa, loc);
    Pose<double, AngleAxis> P1(P);

}

TEST(Pose, PoseInversion)
{
    AngleAxis<double> aa(1, 2, 3);
    double loc[3] = {1, 2, 3};
    Pose<double, AngleAxis> P(aa, loc);
    Pose<double, AngleAxis> inv = P.inverse();
    for (int i = 0; i < 3; ++i)
    {
        ASSERT_DOUBLE_EQ(-aa.data()[i], inv.data()[i]);
        ASSERT_DOUBLE_EQ(-loc[i], inv.data()[i + 3]);
    }
}

TEST(Pose, PoseIdentity)
{
    Pose6d p = Pose6d::Identity();
    Pose34d q = Pose34d::Identity();
}

} // namespace geometry
} // namespace idl

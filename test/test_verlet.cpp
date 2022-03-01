#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class VerletTest : public ::testing::Test
{

protected:
    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->dt = 0.1;

        sys->rx = new double[2]();
        sys->ry = new double[2]();
        sys->rz = new double[2]();
        sys->vx = new double[2]();
        sys->vy = new double[2]();
        sys->vz = new double[2]();
        sys->fx = new double[2]();
        sys->fy = new double[2]();
        sys->fz = new double[2]();
        
        // initial positions 
        sys->rx[0] = -1.0;
        sys->rx[1] = 1.0;

        // initial forces
        sys->fy[0] = 2.0;
        sys->fy[1] = -2.0;
        sys->fz[0] = -2.0;
        sys->fz[1] = 2.0;
    }

    void TearDown()
    {
        delete[] sys->rx;
        delete[] sys->ry;
        delete[] sys->rz;

        delete[] sys->vx;
        delete[] sys->vy;
        delete[] sys->vz;

        delete sys;
    }
};

TEST_F(VerletTest, part1)
{
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
    ASSERT_DOUBLE_EQ(sys->fy[0], 2.0);
    ASSERT_DOUBLE_EQ(sys->fy[1], -2.0);
    ASSERT_DOUBLE_EQ(sys->fz[0], -2.0);
    ASSERT_DOUBLE_EQ(sys->fz[1], 2.0);

    verlet1(sys);

    double exp_v = sys->dt / mvsq2e;
    EXPECT_DOUBLE_EQ(sys->vx[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->vx[1], 0.0);
    EXPECT_DOUBLE_EQ(sys->vy[0], exp_v);
    EXPECT_DOUBLE_EQ(sys->vy[1], -exp_v);
    EXPECT_DOUBLE_EQ(sys->vz[0], -exp_v);
    EXPECT_DOUBLE_EQ(sys->vz[1], exp_v);

    double exp_r = sys->dt * exp_v;
    EXPECT_DOUBLE_EQ(sys->rx[0], -1.0);
    EXPECT_DOUBLE_EQ(sys->rx[1], 1.0);
    EXPECT_DOUBLE_EQ(sys->ry[0], exp_r);
    EXPECT_DOUBLE_EQ(sys->ry[1], -exp_r);
    EXPECT_DOUBLE_EQ(sys->rz[0], -exp_r);
    EXPECT_DOUBLE_EQ(sys->rz[1], exp_r);
}

TEST_F(VerletTest, part2)
{
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
    ASSERT_DOUBLE_EQ(sys->fy[0], 2.0);
    ASSERT_DOUBLE_EQ(sys->fy[1], -2.0);
    ASSERT_DOUBLE_EQ(sys->fz[0], -2.0);
    ASSERT_DOUBLE_EQ(sys->fz[1], 2.0);

    verlet2(sys);

    double exp_v = sys->dt / mvsq2e;
    EXPECT_DOUBLE_EQ(sys->vx[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->vx[1], 0.0);
    EXPECT_DOUBLE_EQ(sys->vy[0], exp_v);
    EXPECT_DOUBLE_EQ(sys->vy[1], -exp_v);
    EXPECT_DOUBLE_EQ(sys->vz[0], -exp_v);
    EXPECT_DOUBLE_EQ(sys->vz[1], exp_v);
}
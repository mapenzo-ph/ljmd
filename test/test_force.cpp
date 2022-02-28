#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class ForceTest: public ::testing::Test 
{

protected:
    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->epsilon = 3.0;
        sys->sigma = 1.0;

        sys->rx = new double[2]();
        sys->ry = new double[2]();
        sys->rz = new double[2]();

        sys->fx = new double[2];
        sys->fy = new double[2];
        sys->fz = new double[2];
        
        // initialize positions
        // take two atoms with dist 2.0 along x
        // this way force is easy to compute
        sys->rx[0] = -1.0;
        sys->rx[1] = 1.0;
    }

    void TearDown()
    {
        delete[] sys->rx;
        delete[] sys->ry;
        delete[] sys->rz;

        delete[] sys->fx;
        delete[] sys->fy;
        delete[] sys->fz;

        delete sys;
    }
};

TEST_F(ForceTest, inside_rcut_direct)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
    ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
    ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

    sys->box = 10.0;
    sys->rcut = 3.0;
    force(sys);

    double exp_f = -sys->epsilon * 93.0 / 512.0;
    EXPECT_DOUBLE_EQ(sys->fx[0], -exp_f);
    EXPECT_DOUBLE_EQ(sys->fx[1], exp_f);
    EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
    EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);

    double exp_epot = -sys->epsilon * 63.0 / 1024.0;
    EXPECT_DOUBLE_EQ(sys->epot, exp_epot);
}

TEST_F(ForceTest, inside_rcut_pbc)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
    ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
    ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

    sys->box = 3.0;
    sys->rcut = 1.5;
    force(sys);

    double exp_f = 24.0 * sys->epsilon;
    EXPECT_DOUBLE_EQ(sys->fx[0], exp_f);
    EXPECT_DOUBLE_EQ(sys->fx[1], -exp_f);
    EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
    EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);

    EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}

TEST_F(ForceTest, outside_rcut_direct)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
    ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
    ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

    sys->box = 10.0;
    sys->rcut = 1.0;
    force(sys);

    EXPECT_DOUBLE_EQ(sys->fx[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fx[1], 0.0);
    EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
    EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
    EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);
    
    EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}



#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class ForceTest: public ::testing::Test 
{

protected:
    mdsys_t *sys;
    particles_t *part;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->epsilon = 3.0;
        sys->sigma = 1.0;

        part = new particles_t[2]();
        
        // initialize positions
        // take two atoms with dist 2.0 along x
        // this way force is easy to compute
        part[0].rx = -1.0;
        part[1].rx = 1.0;
    }

    void TearDown()
    {
        delete[] part;

        delete sys;
    }
};

TEST_F(ForceTest, inside_rcut_direct)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(part[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(part[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(part[0].ry, 0.0);
    ASSERT_DOUBLE_EQ(part[1].ry, 0.0);
    ASSERT_DOUBLE_EQ(part[0].rz, 0.0);
    ASSERT_DOUBLE_EQ(part[1].rz, 0.0);

    sys->box = 10.0;
    sys->rcut = 3.0;
    force(sys,part);

    double exp_f = -sys->epsilon * 93.0 / 512.0;
    EXPECT_DOUBLE_EQ(part[0].fx, -exp_f);
    EXPECT_DOUBLE_EQ(part[1].fx, exp_f);
    EXPECT_DOUBLE_EQ(part[0].fy, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fy, 0.0);
    EXPECT_DOUBLE_EQ(part[0].fz, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fz, 0.0);

    double exp_epot = -sys->epsilon * 63.0 / 1024.0;
    EXPECT_DOUBLE_EQ(sys->epot, exp_epot);
}

TEST_F(ForceTest, inside_rcut_pbc)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(part[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(part[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(part[0].ry, 0.0);
    ASSERT_DOUBLE_EQ(part[1].ry, 0.0);
    ASSERT_DOUBLE_EQ(part[0].rz, 0.0);
    ASSERT_DOUBLE_EQ(part[1].rz, 0.0);

    sys->box = 3.0;
    sys->rcut = 1.5;
    force(sys,part);

    double exp_f = 24.0 * sys->epsilon;
    EXPECT_DOUBLE_EQ(part[0].fx, exp_f);
    EXPECT_DOUBLE_EQ(part[1].fx, -exp_f);
    EXPECT_DOUBLE_EQ(part[0].fy, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fy, 0.0);
    EXPECT_DOUBLE_EQ(part[0].fz, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fz, 0.0);

    EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}

TEST_F(ForceTest, outside_rcut_direct)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(part[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(part[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(part[0].ry, 0.0);
    ASSERT_DOUBLE_EQ(part[1].ry, 0.0);
    ASSERT_DOUBLE_EQ(part[0].rz, 0.0);
    ASSERT_DOUBLE_EQ(part[1].rz, 0.0);

    sys->box = 10.0;
    sys->rcut = 1.0;
    force(sys,part);

    EXPECT_DOUBLE_EQ(part[0].fx, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fx, 0.0);
    EXPECT_DOUBLE_EQ(part[0].fy, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fy, 0.0);
    EXPECT_DOUBLE_EQ(part[0].fz, 0.0);
    EXPECT_DOUBLE_EQ(part[1].fz, 0.0);
    
    EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}



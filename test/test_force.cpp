#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}


class ForceTest: public ::testing::Test 
{

public:
    mdsys_t *sys;
    coords_t *coord;
    for_t *forces;

    void SetUp()
    {

        sys = new mdsys_t;

        sys->natoms = 2;
        sys->epsilon = 3.0;
        sys->sigma = 1.0;

        coord = new coords_t[2]();
        forces = new for_t[2]();
        
        // initialize positions
        // take two atoms with dist 2.0 along x
        // this way force is easy to compute
        coord[0].rx = -1.0;
        coord[1].rx = 1.0;
    }

    void TearDown()
    {
        delete[] coord;
        delete[] forces;
        delete sys;

    }
};

TEST_F(ForceTest, inside_rcut_direct)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(coord[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(coord[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(coord[0].ry, 0.0);
    ASSERT_DOUBLE_EQ(coord[1].ry, 0.0);
    ASSERT_DOUBLE_EQ(coord[0].rz, 0.0);
    ASSERT_DOUBLE_EQ(coord[1].rz, 0.0);

    sys->box = 10.0;
    sys->rcut = 3.0;

    force(sys,coord,forces);

    double exp_f = -sys->epsilon * 93.0 / 512.0;
    EXPECT_DOUBLE_EQ(forces[0].fx, -exp_f);
    EXPECT_DOUBLE_EQ(forces[1].fx, exp_f);
    EXPECT_DOUBLE_EQ(forces[0].fy, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fy, 0.0);
    EXPECT_DOUBLE_EQ(forces[0].fz, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fz, 0.0);

    double exp_epot = -sys->epsilon * 63.0 / 1024.0;
    EXPECT_DOUBLE_EQ(sys->epot, exp_epot);
}

TEST_F(ForceTest, inside_rcut_pbc)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(coord[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(coord[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(coord[0].ry, 0.0);
    ASSERT_DOUBLE_EQ(coord[1].ry, 0.0);
    ASSERT_DOUBLE_EQ(coord[0].rz, 0.0);
    ASSERT_DOUBLE_EQ(coord[1].rz, 0.0);

    sys->box = 3.0;
    sys->rcut = 1.5;

    force(sys,coord,forces);

    double exp_f = 24.0 * sys->epsilon;
    EXPECT_DOUBLE_EQ(forces[0].fx, exp_f);
    EXPECT_DOUBLE_EQ(forces[1].fx, -exp_f);
    EXPECT_DOUBLE_EQ(forces[0].fy, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fy, 0.0);
    EXPECT_DOUBLE_EQ(forces[0].fz, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fz, 0.0);

    EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}

TEST_F(ForceTest, outside_rcut_direct)
{
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(coord[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(coord[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(coord[0].ry, 0.0);
    ASSERT_DOUBLE_EQ(coord[1].ry, 0.0);
    ASSERT_DOUBLE_EQ(coord[0].rz, 0.0);
    ASSERT_DOUBLE_EQ(coord[1].rz, 0.0);

    sys->box = 10.0;
    sys->rcut = 1.0;

    force(sys,coord,forces);

    EXPECT_DOUBLE_EQ(forces[0].fx, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fx, 0.0);
    EXPECT_DOUBLE_EQ(forces[0].fy, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fy, 0.0);
    EXPECT_DOUBLE_EQ(forces[0].fz, 0.0);
    EXPECT_DOUBLE_EQ(forces[1].fz, 0.0);
    
    EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}





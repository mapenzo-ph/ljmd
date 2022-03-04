#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class VerletTest : public ::testing::Test
{

protected:
    mdsys_t *sys;
    coords_t *coord;
    for_t *forces;
    vel_t *vel;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->dt = 0.1;
        
        coord = new coords_t[2]();
        forces = new for_t[2]();
        vel = new vel_t[2]();
        
        
        // initial positions 
        coord[0].rx = -1.0;
        coord[1].rx = 1.0;

        // initial forces
        forces[0].fy = 2.0;
        forces[1].fy = -2.0;
        forces[0].fz = -2.0;
        forces[1].fz = 2.0;
    }

    void TearDown()
    {
        delete[] coord;
        delete[] forces;
        delete[] vel;

        delete sys;
    }
};

TEST_F(VerletTest, part1)
{
    ASSERT_DOUBLE_EQ(coord[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(coord[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(forces[0].fy, 2.0);
    ASSERT_DOUBLE_EQ(forces[1].fy, -2.0);
    ASSERT_DOUBLE_EQ(forces[0].fz, -2.0);
    ASSERT_DOUBLE_EQ(forces[1].fz, 2.0);

    verlet1(sys,coord,vel,forces);

    double exp_v = sys->dt / mvsq2e;
    EXPECT_DOUBLE_EQ(vel[0].vx, 0.0);
    EXPECT_DOUBLE_EQ(vel[1].vx, 0.0);
    EXPECT_DOUBLE_EQ(vel[0].vy, exp_v);
    EXPECT_DOUBLE_EQ(vel[1].vy, -exp_v);
    EXPECT_DOUBLE_EQ(vel[0].vz, -exp_v);
    EXPECT_DOUBLE_EQ(vel[1].vz, exp_v);

    double exp_r = sys->dt * exp_v;
    EXPECT_DOUBLE_EQ(coord[0].rx, -1.0);
    EXPECT_DOUBLE_EQ(coord[1].rx, 1.0);
    EXPECT_DOUBLE_EQ(coord[0].ry, exp_r);
    EXPECT_DOUBLE_EQ(coord[1].ry, -exp_r);
    EXPECT_DOUBLE_EQ(coord[0].rz, -exp_r);
    EXPECT_DOUBLE_EQ(coord[1].rz, exp_r);
}

TEST_F(VerletTest, part2)
{
    ASSERT_DOUBLE_EQ(coord[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(coord[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(forces[0].fy, 2.0);
    ASSERT_DOUBLE_EQ(forces[1].fy, -2.0);
    ASSERT_DOUBLE_EQ(forces[0].fz, -2.0);
    ASSERT_DOUBLE_EQ(forces[1].fz, 2.0);

    verlet2(sys,vel,forces);

    double exp_v = sys->dt / mvsq2e;
    EXPECT_DOUBLE_EQ(vel[0].vx, 0.0);
    EXPECT_DOUBLE_EQ(vel[1].vx, 0.0);
    EXPECT_DOUBLE_EQ(vel[0].vy, exp_v);
    EXPECT_DOUBLE_EQ(vel[1].vy, -exp_v);
    EXPECT_DOUBLE_EQ(vel[0].vz, -exp_v);
    EXPECT_DOUBLE_EQ(vel[1].vz, exp_v);
}
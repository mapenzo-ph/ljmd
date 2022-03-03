#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class VerletTest : public ::testing::Test
{

protected:
    mdsys_t *sys;
    particles_t *part;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->dt = 0.1;
        
        part = new particles_t[2]();
        
        // initial positions 
        part[0].rx = -1.0;
        part[1].rx = 1.0;

        // initial forces
        part[0].fy = 2.0;
        part[1].fy = -2.0;
        part[0].fz = -2.0;
        part[1].fz = 2.0;
    }

    void TearDown()
    {
        delete[] part;

        delete sys;
    }
};

TEST_F(VerletTest, part1)
{
    ASSERT_DOUBLE_EQ(part[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(part[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(part[0].fy, 2.0);
    ASSERT_DOUBLE_EQ(part[1].fy, -2.0);
    ASSERT_DOUBLE_EQ(part[0].fz, -2.0);
    ASSERT_DOUBLE_EQ(part[1].fz, 2.0);

    verlet1(sys,part);

    double exp_v = sys->dt / mvsq2e;
    EXPECT_DOUBLE_EQ(part[0].vx, 0.0);
    EXPECT_DOUBLE_EQ(part[1].vx, 0.0);
    EXPECT_DOUBLE_EQ(part[0].vy, exp_v);
    EXPECT_DOUBLE_EQ(part[1].vy, -exp_v);
    EXPECT_DOUBLE_EQ(part[0].vz, -exp_v);
    EXPECT_DOUBLE_EQ(part[1].vz, exp_v);

    double exp_r = sys->dt * exp_v;
    EXPECT_DOUBLE_EQ(part[0].rx, -1.0);
    EXPECT_DOUBLE_EQ(part[1].rx, 1.0);
    EXPECT_DOUBLE_EQ(part[0].ry, exp_r);
    EXPECT_DOUBLE_EQ(part[1].ry, -exp_r);
    EXPECT_DOUBLE_EQ(part[0].rz, -exp_r);
    EXPECT_DOUBLE_EQ(part[1].rz, exp_r);
}

TEST_F(VerletTest, part2)
{
    ASSERT_DOUBLE_EQ(part[0].rx, -1.0);
    ASSERT_DOUBLE_EQ(part[1].rx, 1.0);
    ASSERT_DOUBLE_EQ(part[0].fy, 2.0);
    ASSERT_DOUBLE_EQ(part[1].fy, -2.0);
    ASSERT_DOUBLE_EQ(part[0].fz, -2.0);
    ASSERT_DOUBLE_EQ(part[1].fz, 2.0);

    verlet2(sys,part);

    double exp_v = sys->dt / mvsq2e;
    EXPECT_DOUBLE_EQ(part[0].vx, 0.0);
    EXPECT_DOUBLE_EQ(part[1].vx, 0.0);
    EXPECT_DOUBLE_EQ(part[0].vy, exp_v);
    EXPECT_DOUBLE_EQ(part[1].vy, -exp_v);
    EXPECT_DOUBLE_EQ(part[0].vz, -exp_v);
    EXPECT_DOUBLE_EQ(part[1].vz, exp_v);
}
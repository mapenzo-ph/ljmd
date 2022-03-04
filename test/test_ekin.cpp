#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class EKinTest : public ::testing::Test
{

protected:
    mdsys_t *sys;
    vel_t *vel;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 3;
        sys->mass = 1.0;

        vel = new vel_t[3]();

        // take velocities along axes for simplicity
        vel[0].vx = 1.0;
        vel[1].vy = 3.0;
        vel[2].vz = -2.0;
    }

    void TearDown()
    {
        delete[] vel;

        delete sys;
    }
};

TEST_F(EKinTest, basic_test)
{
    ASSERT_DOUBLE_EQ(vel[0].vx, 1.0);
    ASSERT_DOUBLE_EQ(vel[1].vx, 0.0);
    ASSERT_DOUBLE_EQ(vel[2].vx, 0.0);
    ASSERT_DOUBLE_EQ(vel[0].vy, 0.0);
    ASSERT_DOUBLE_EQ(vel[1].vy, 3.0);
    ASSERT_DOUBLE_EQ(vel[2].vy, 0.0);
    ASSERT_DOUBLE_EQ(vel[0].vz, 0.0);
    ASSERT_DOUBLE_EQ(vel[1].vz, 0.0);
    ASSERT_DOUBLE_EQ(vel[2].vz, -2.0);

    ekin(sys,vel);

    double exp_ekin = 7.0*mvsq2e;
    EXPECT_DOUBLE_EQ(sys->ekin, exp_ekin);

    double exp_temp = exp_ekin / 3.0 / kboltz;
    EXPECT_DOUBLE_EQ(sys->temp, exp_temp);
}

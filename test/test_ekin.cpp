#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class EKinTest : public ::testing::Test
{

protected:
    mdsys_t *sys;
    particles_t *part;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 3;
        sys->mass = 1.0;

        part = new particles_t[3]();

        // take velocities along axes for simplicity
        part[0].vx = 1.0;
        part[1].vy = 3.0;
        part[2].vz = -2.0;
    }

    void TearDown()
    {
        delete[] part;

        delete sys;
    }
};

TEST_F(EKinTest, basic_test)
{
    ASSERT_DOUBLE_EQ(part[0].vx, 1.0);
    ASSERT_DOUBLE_EQ(part[1].vx, 0.0);
    ASSERT_DOUBLE_EQ(part[2].vx, 0.0);
    ASSERT_DOUBLE_EQ(part[0].vy, 0.0);
    ASSERT_DOUBLE_EQ(part[1].vy, 3.0);
    ASSERT_DOUBLE_EQ(part[2].vy, 0.0);
    ASSERT_DOUBLE_EQ(part[0].vz, 0.0);
    ASSERT_DOUBLE_EQ(part[1].vz, 0.0);
    ASSERT_DOUBLE_EQ(part[2].vz, -2.0);

    ekin(sys,part);

    double exp_ekin = 7.0*mvsq2e;
    EXPECT_DOUBLE_EQ(sys->ekin, exp_ekin);

    double exp_temp = exp_ekin / 3.0 / kboltz;
    EXPECT_DOUBLE_EQ(sys->temp, exp_temp);
}

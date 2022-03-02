#include "gtest/gtest.h"
extern "C" {
#include "utils.h"
}

class EKinTest : public ::testing::Test
{

protected:
    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 3;
        sys->mass = 1.0;

        sys->vx = new double[3]();
        sys->vy = new double[3]();
        sys->vz = new double[3]();

        // take velocities along axes for simplicity
        sys->vx[0] = 1.0;
        sys->vy[1] = 3.0;
        sys->vz[2] = -2.0;
    }

    void TearDown()
    {
        delete[] sys->vx;
        delete[] sys->vy;
        delete[] sys->vz;

        delete sys;
    }
};

TEST_F(EKinTest, basic_test)
{
    ASSERT_DOUBLE_EQ(sys->vx[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->vx[1], 0.0);
    ASSERT_DOUBLE_EQ(sys->vx[2], 0.0);
    ASSERT_DOUBLE_EQ(sys->vy[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->vy[1], 3.0);
    ASSERT_DOUBLE_EQ(sys->vy[2], 0.0);
    ASSERT_DOUBLE_EQ(sys->vz[0], 0.0);
    ASSERT_DOUBLE_EQ(sys->vz[1], 0.0);
    ASSERT_DOUBLE_EQ(sys->vz[2], -2.0);

    ekin(sys);

    double exp_ekin = 7.0*mvsq2e;
    EXPECT_DOUBLE_EQ(sys->ekin, exp_ekin);

    double exp_temp = exp_ekin / 3.0 / kboltz;
    EXPECT_DOUBLE_EQ(sys->temp, exp_temp);
}

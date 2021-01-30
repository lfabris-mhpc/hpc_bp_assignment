// -*- c++ -*-

#include <gtest/gtest.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

extern "C" {
#include <utilities.h>
}

TEST(UtilitiesWallclock, elapsed)
{
	const double wait = 0.1;
	
	struct timespec req, rem;
	req.tv_sec = 0;
	req.tv_nsec = wait * 1000 * 1000 * 1000;
	
	double delta = -wallclock();
	nanosleep(&req, &rem);
	delta += wallclock();
	
	ASSERT_GE(delta, wait);
	ASSERT_NEAR(delta, wait, 0.001);
}

TEST(UtilitiesAzzero, initAndAzzero)
{
	constexpr int len{32};
	double v[len];
	for (auto i = 0; i < len; ++i) {
		v[i] = (double) (i + 1);
	}

	for (auto i = 0; i < len; ++i) {
		ASSERT_DOUBLE_EQ(v[i], (double) (i + 1));
	}

	azzero(v, len);

	for (auto i = 0; i < len; ++i) {
		ASSERT_DOUBLE_EQ(v[i], 0.0);
	}
}

TEST(UtilitiesPBC, sequentialPBC)
{
	constexpr double radius_inner{1.1};
	for (auto i = -10; i < 11; ++i) {
		double toss{0.2 * i};
		double contracted{pbc(toss, radius_inner)};
	
		ASSERT_GE(contracted, -radius_inner);
		ASSERT_LE(contracted, radius_inner);
	}
}

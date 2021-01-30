// -*- c++ -*-

#include <gtest/gtest.h>
#include <time.h>
#include <sys/time.h>

extern "C" {
#include <utilities.h>
}

TEST(Wallclock, elapsed)
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
/*
TEST(Color, presets)
{
	ASSERT_EQ(white.r,255);
	ASSERT_EQ(white.g,255);
	ASSERT_EQ(white.b,255);
	ASSERT_EQ(black.r,0);
	ASSERT_EQ(black.g,0);
	ASSERT_EQ(black.b,0);
	ASSERT_EQ(red.r,255);
	ASSERT_EQ(red.g,0);
	ASSERT_EQ(red.b,0);
	ASSERT_EQ(green.r,0);
	ASSERT_EQ(green.g,255);
	ASSERT_EQ(green.b,0);
	ASSERT_EQ(blue.r,0);
	ASSERT_EQ(blue.g,0);
	ASSERT_EQ(blue.b,255);
	ASSERT_EQ(yellow.r,255);
	ASSERT_EQ(yellow.g,255);
	ASSERT_EQ(yellow.b,0);
}
*/

#pragma once

#include <cmath>
#include <algorithm>
#include <random>

namespace astrokit
{

	inline double safe_acos(double x)
	//clamps to -1.0 <= x <= 1.0 then computes the acos
	{
		x = std::clamp(x, -1.0, 1.0); //make sure there's no machine precision/roundoff error extending beyond the acos bounds
		return std::acos(x);
	}

	inline int random_int(int min, int max)
	//using a uniform distribution for these
	{
		static std::mt19937 rng(std::random_device{}());   // only seeded once
		std::uniform_int_distribution<int> dist(min, max);
		return dist(rng);
	}

	inline double random_double(double min, double max)
	{
		static std::mt19937 rng(std::random_device{}());
		std::uniform_real_distribution<double> dist(min, max);
		return dist(rng);
	}

} // namespace astrokit



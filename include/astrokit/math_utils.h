#pragma once

#include <cmath>
#include <algorithm>
#include <random>
#include <Eigen/Dense>
#include "constants.h"

namespace astrokit
{

	inline double safe_acos(double x)
	//clamps to -1.0 <= x <= 1.0 then computes the acos
	{
		x = std::clamp(x, -1.0, 1.0); //make sure there's no machine precision/roundoff error extending beyond the acos bounds
		return std::acos(x);
	}

	inline int sign(double x)
	{
		//useful utility function that returns:
		//            / +1 if x > 0
		// sign(x) = -   0 if x = 0
		//            \ -1 if x < 0
		return (x > 0) - (x < 0);
	}

	inline double angle_between_vecs(Eigen::Vector3d vec1, Eigen::Vector3d vec2)
	//directionless; simple angle between two vectors
	{
		return safe_acos(vec1.dot(vec2) / (vec1.norm() * vec2.norm()));
	}

	inline double angle_between_vecs_w_direction(Eigen::Vector3d vec1, Eigen::Vector3d vec2, Eigen::Vector3d ref_axis)
	//uses the reference axis to determine if the rotation angle is > or < 180 degrees
	{
		//first get the angle
		double angle = angle_between_vecs(vec1, vec2);
		
		//then check vec1 x vec2 against ref_axis
		Eigen::Vector3d xprod = vec1.cross(vec2);
		double xprod_angle = angle_between_vecs(ref_axis, xprod);
		if (xprod_angle > astrokit::PI / 2.0)
		{
			return 2.0 * astrokit::PI - angle;
		}
		return angle;
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



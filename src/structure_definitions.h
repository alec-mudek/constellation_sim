
#pragma once
#include <Eigen/Dense>

struct State //Contains both the cartesian (ICRF) state and the orbital elements for each time step, et
{
	double et;
	Eigen::Vector3d pos;
	Eigen::Vector3d vel;
	double sma;
	double ecc;
	double inc;
	double raan;
	double argp;
	double ta;
};

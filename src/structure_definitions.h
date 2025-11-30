
#pragma once
#include <Eigen/Dense>

struct State //Contains both the cartesian (ICRF) state and the orbital elements for each time step, et
{
	double et;

	//icrf cartesian state
	Eigen::Vector3d pos;
	Eigen::Vector3d vel;
	
	//coes
	double sma;
	double ecc;
	double inc;
	double raan;
	double argp;
	double ta;
};

struct TrackingState
{
	double et;

	//need to know relevant mean orbital elements for bounding box
	double sma_mean;
	double ecc_mean;
	double inc_mean;
	double raan_mean;
	double arglat_mean; //instead of argp + ta

	//need to know relative position in constellation
	double neighbor_one_rel_ta; //true anomaly delta b/w a spacecraft and its neighbor
	double neighbor_two_rel_ta;
};

struct COE //using this to store reference conic trajectories
{
	double sma;
	double ecc;
	double inc;
	double raan;
	double argp;
	double ta;
};

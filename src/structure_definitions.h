
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
	double inc_mean;
	double raan_mean;

	//need to know relative position in constellation
	double neighbor1_rel_angle; //approximate argument of latitude delta b/w a spacecraft and its neighbor (actually find the angle between the position vectors to avoid angle wrap headaches)
	double neighbor2_rel_angle;
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

struct BoundingBox //orbital element bounding boxes for members of constellation
//note: using orbital element bounds instead of RTN for this simulation
{
	double dsma;
	double dinc;
	double draan;
	double darglat;
};

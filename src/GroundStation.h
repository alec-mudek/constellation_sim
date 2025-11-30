
#pragma once
#include <string>
#include <astrokit/constants.h>
#include <Eigen/Dense>

class GroundStation
{
public:
	GroundStation(std::string name);
	GroundStation(std::string name, double lon, double lat);
	GroundStation(std::string name, double lon, double lat, double elev_mask);
	
	//getters
	std::string get_name() const;
	double get_lon() const;
	double get_lat() const;
	double get_elevation_mask() const;

	Eigen::Vector3d get_surface_normal_bcf() const;

	//setters
	void set_name(std::string new_name);
	void set_lon(double new_lon);
	void set_lat(double new_lat);
	void set_elevation_mask(double new_elev_mask);

	//utilities
	void compute_surface_normal();

private:
	std::string name;

	double lon;
	double lat;

	Eigen::Vector3d bcf_pos; //cartesian bcf position vector
	Eigen::Vector3d bcf_surf_normal; //unit surface normal vector

	double elevation_mask; //elevation angle at which the spacraft can be considered "in-range" of ground station (must be >= mask)
	
};


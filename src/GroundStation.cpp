#include "GroundStation.h"


GroundStation::GroundStation(std::string name) :
	lon(0.0), lat(0.0), elevation_mask(0.0)
{
	set_name(name);
}

GroundStation::GroundStation(std::string name, double lon, double lat) : elevation_mask(0.0)
{
	set_name(name);
	set_lon(lon);
	set_lat(lat);
}

GroundStation::GroundStation(std::string name, double lon, double lat, double elev_mask)
{
	set_name(name);
	set_lon(lon);
	set_lat(lat);
	set_elevation_mask(elev_mask);
}

#pragma region getters
std::string GroundStation::get_name() const
{
	return this->name;
}
double GroundStation::get_lon() const
{
	return this->lon;
}

double GroundStation::get_lat() const
{
	return this->lat;
}

double GroundStation::get_elevation_mask() const
{
	return this->elevation_mask;
}
Eigen::Vector3d GroundStation::get_surface_normal_bcf() const
{
	return this->bcf_surf_normal;
}
#pragma endregion getters

#pragma region setters
void GroundStation::set_name(std::string new_name)
{
	this->name = new_name;
}
void GroundStation::set_lon(double new_lon)
{
	this->lon = new_lon;
}

void GroundStation::set_lat(double new_lat)
{
	this->lat = new_lat;
}

void GroundStation::set_elevation_mask(double new_elev_mask)
{
	this->elevation_mask = new_elev_mask;
}
#pragma endregion setters

#pragma region utilities
void GroundStation::compute_surface_normal()
{
	//currently modeling a perfect sphere for the planet (which is why this logic can live in the ground station class)
	//the surface normal is just the unit of the ground station position vector which we know from the lat/lon
	this->bcf_surf_normal = Eigen::Vector3d(std::cos(lon) * std::cos(lat), std::sin(lon) * std::cos(lat), std::sin(lat));
	this->bcf_surf_normal.normalize(); //just to be safe
}
#pragma endregion utilities

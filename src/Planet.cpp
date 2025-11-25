#include "Planet.h"

Planet::Planet(SpiceHandler& spice) :
	mu(1.0), mean_radius(1.0), eq_radius(1.0), j2(0.0), spkid(0), bcf_frame_name(""), stations{}, spice(spice)
{
}

Planet::Planet(SpiceHandler& spice, double mu, double mean_radius, double eq_radius, double j2, int spkid, std::string bcf_frame_name) :
	spice(spice), stations{}
{
	set_mu(mu);
	set_mean_radius(mean_radius);
	set_eq_radius(eq_radius);
	set_j2(j2);
	set_spkid(spkid);
	set_bcf_frame_name(bcf_frame_name);
}



Planet::Planet(SpiceHandler& spice, double mu, double mean_radius, double eq_radius, double j2, int spkid, std::string bcf_frame_name, std::vector<std::array<double, 2>> station_lon_lats) :
	spice(spice)
{
	set_mu(mu);
	set_mean_radius(mean_radius);
	set_eq_radius(eq_radius);
	set_j2(j2);
	set_spkid(spkid);
	set_bcf_frame_name(bcf_frame_name);
	for (std::array<double, 2> lon_lat : station_lon_lats)
	{
		new_station(lon_lat[0], lon_lat[1]);
		//note: this constructor doesn't allow for specifying the half_angle of the stations
		//      if you want that level of detail you need to use the correct NewStation() method directly to add them
		//      ground stations default to a 90 degree half angle
	}
}

#pragma region getters
double Planet::get_mu() const
{
	return this->mu;
}

double Planet::get_mean_radius() const
{
	return this->mean_radius;
}

double Planet::get_eq_radius() const
{
	return this->eq_radius;
}

double Planet::get_j2() const
{
	return this->j2;
}

int Planet::get_spkid() const
{
	return this->spkid;
}

std::string Planet::get_bcf_frame_name() const
{
	return this->bcf_frame_name;
}

const std::vector<GroundStation>& Planet::get_stations() const
{
	return this->stations;
}
#pragma endregion getters

#pragma region setters
void Planet::set_mu(double new_mu)
{
	this->mu = new_mu;
}

void Planet::set_mean_radius(double new_radius)
{
	this->mean_radius = new_radius;
}

void Planet::set_eq_radius(double new_radius)
{
	this->eq_radius = new_radius;
}

void Planet::set_j2(double new_j2)
{
	this->j2 = new_j2;
}

void Planet::set_spkid(int new_spkid)
{
	this->spkid = new_spkid;
}

void Planet::set_bcf_frame_name(std::string new_frame_name)
{
	this->bcf_frame_name = new_frame_name;
}
#pragma endregion setters

#pragma region utilities
void Planet::new_station(double lon, double lat)
{
	GroundStation new_station(lon, lat);
	this->stations.push_back(new_station);
}

void Planet::new_station(double lon, double lat, double half_angle)
{
	GroundStation new_station(lon, lat, half_angle);
	this->stations.push_back(new_station);
}

void Planet::remove_station(double lon, double lat)
{
	for (int i = 0; i < this->stations.size(); i++)
	{
		GroundStation gs = this->stations[i];
		if (gs.get_lon() == lon && gs.get_lat() == lat)
		{
			stations.erase(stations.begin() + i); //assumes only 1 station per lat/lon; will only remove the first if there are duplicates
			break;
		}
	}
}

Eigen::Vector3d Planet::sun_vector(double et)
{
	Eigen::Vector3d r_sun = spice.fetch_pos(et, 10, this->spkid, "J2000");

	return r_sun.normalized();
}

Eigen::Matrix3d Planet::icrf_R_bcf(double et)
{
	Eigen::Matrix3d icrf_C_bcf = spice.fetch_rot_matrix(et, "J2000", this->bcf_frame_name);
	
	return icrf_C_bcf;
}

Eigen::Matrix3d Planet::bcf_R_icrf(double et)
{
	Eigen::Matrix3d bcf_C_icrf = icrf_R_bcf(et).transpose();
	
	return bcf_C_icrf;
}

#pragma endregion utilities

#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "SpiceHandler.h"
#include "GroundStation.h"

class Planet
{
public:
	Planet(SpiceHandler& spice);
	Planet(SpiceHandler& spice, double mu, double mean_radius, double eq_radius, double j2, int spkid, std::string bcf_frame_name);
	Planet(SpiceHandler& spice, double mu, double mean_radius, double eq_radius, double j2, int spkid, std::string bcf_frame_name, 
		   std::vector<std::string> station_names, std::vector<std::array<double, 2>> station_lon_lats);

	//going for a singleton-ish pattern for the Planet class; don't want it to be copyable
	Planet(const Planet&) = delete;
	Planet& operator=(const Planet&) = delete;
	Planet(Planet&&) = delete;
	Planet& operator=(Planet&&) = delete;

	//getters
	double get_mu() const;
	double get_mean_radius() const;
	double get_eq_radius() const;
	double get_j2() const;
	int get_spkid() const;
	std::string get_bcf_frame_name() const;
	const std::vector<GroundStation>& get_stations() const;

	//setters
	void set_mu(double new_mu);
	void set_mean_radius(double new_radius);
	void set_eq_radius(double new_radius);
	void set_j2(double new_j2);
	void set_spkid(int new_spkid);
	void set_bcf_frame_name(std::string new_frame_name);

	//utilities
	void new_station(std::string name, double lon, double lat);
	void new_station(std::string name, double lon, double lat, double half_angle);
	void remove_station(double lon, double lat); //assumes only 1 station per lat/lon location
	
	Eigen::Vector3d sun_vector(double et); //unit vector from the center of the planet to the sun
	//note: rotation matrices follow the convection v_new = C * v_old
	Eigen::Matrix3d icrf_R_bcf(double et);
	Eigen::Matrix3d bcf_R_icrf(double et);
	
private:
	double mu;
	double mean_radius; //assuming a perfect spheroid for the body model (lat/lon computations & altitude)
	double j2;     //assuming an oblate spheroid for the gravity modeling
	double eq_radius; //need the equatorial radius for J2 EOMs
	//note: future work; model an oblate spheroid or upgrade to a full obj file & higher-order gravity

	int spkid; //spice id
	std::string bcf_frame_name; //e.g. 'IAU_EARTH', 'IAU_MARS', etc.; the spice name for the planet's BCF frame
	//note: there are more precise frames than the IAU, especially for earth; just using the IAU for this demo

	std::vector<GroundStation> stations; //variable number of ground stations

	SpiceHandler& spice; //just a reference here; singleton pattern for SpiceHandler
};


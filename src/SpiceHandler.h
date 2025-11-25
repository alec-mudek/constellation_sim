
#pragma once
#include <string>
#include <cspice/SpiceUsr.h>
#include <Eigen/Dense>


class SpiceHandler
{
public: 
	SpiceHandler();
	SpiceHandler(std::string de_path, std::string naif_path, std::string pck_path);

	//going for a singleton-ish pattern for the SpiceHandler class; don't want it to be copyable
	SpiceHandler(const SpiceHandler&) = delete;
	SpiceHandler& operator=(const SpiceHandler&) = delete;
	SpiceHandler(SpiceHandler&&) = delete;
	SpiceHandler& operator=(SpiceHandler&&) = delete;

	//getters
	std::string get_de_path() const;
	std::string get_naif_path() const;
	std::string get_pck_path() const;

	//setters
	void set_de_path(std::string new_de_path);
	void set_naif_path(std::string new_naif_path);
	void set_pck_path(std::string new_pck_path);

	//utilities
	void load_kernels() const;
	void reload_kernels() const;
	Eigen::Vector3d fetch_pos(double et, int target_spkid, int observer_spkid, std::string frame_name); 
	//              ^ retrieve cartesian position from spice data
	std::array<Eigen::Vector3d, 2> fetch_state(double et, int target_spkid, int observer_spkid, std::string frame_name); 
	//                             ^ retrieve full carteisan state from spice data
	Eigen::Matrix3d fetch_rot_matrix(double et, std::string from_frame, std::string to_frame);
	double str_date_to_et(std::string date_string);


private:
	std::string de_path; //path to de4XX.bsp file (defaults to de440s.bsp)
	std::string naif_path; //path to naif0012.tls file
	std::string pck_path; //path to pck00011.tpc file

};


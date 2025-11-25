#include "SpiceHandler.h"


SpiceHandler::SpiceHandler() :
	de_path("../kernels/de440s.bsp"),
	naif_path("../kernels/naif0012.tls"),
	pck_path("../kernels/pck00011.tpc")
{
	load_kernels();
}

SpiceHandler::SpiceHandler(std::string de_path, std::string naif_path, std::string pck_path)
{
	set_de_path(de_path);
	set_naif_path(naif_path);
	set_pck_path(pck_path);

	load_kernels();
}

#pragma region getters
std::string SpiceHandler::get_de_path() const
{
	return this->de_path;
}

std::string SpiceHandler::get_naif_path() const
{
	return this->naif_path;
}

std::string SpiceHandler::get_pck_path() const
{
	return this->pck_path;
}
#pragma endregion getters

#pragma region setters
void SpiceHandler::set_de_path(std::string new_de_path)
{
	this->de_path = new_de_path;
}

void SpiceHandler::set_naif_path(std::string new_naif_path)
{
	this->naif_path = new_naif_path;
}

void SpiceHandler::set_pck_path(std::string new_pck_path)
{
	this->pck_path = new_pck_path;
}
#pragma endregion setters

#pragma region utilities
void SpiceHandler::load_kernels() const
{
	furnsh_c(this->de_path.c_str());
	furnsh_c(this->naif_path.c_str());
	furnsh_c(this->pck_path.c_str());
}

void SpiceHandler::reload_kernels() const
{
	kclear_c();
	load_kernels();
}

Eigen::Vector3d SpiceHandler::fetch_pos(double et, int target_spkid, int observer_spkid, std::string frame_name)
{
	std::array<Eigen::Vector3d, 2> state = fetch_state(et, target_spkid, observer_spkid, frame_name);
	return state[0];
}

std::array<Eigen::Vector3d, 2> SpiceHandler::fetch_state(double et, int target_spkid, int observer_spkid, std::string frame_name)
{
	SpiceDouble spice_state[6];
	SpiceDouble lt;
	spkez_c(target_spkid, et, frame_name.c_str(), "None", observer_spkid, spice_state, &lt);
	
	Eigen::Vector3d pos{ spice_state[0], spice_state[1], spice_state[2] };
	Eigen::Vector3d vel{ spice_state[3], spice_state[4], spice_state[5] };

	return std::array<Eigen::Vector3d, 2>{pos, vel};
}
Eigen::Matrix3d SpiceHandler::fetch_rot_matrix(double et, std::string from_frame, std::string to_frame)
{
	//pxform requires a SpiceDouble object for the rotation matrix, C
	//create one to retrieve the rot matrix information, then transfer the data
	// to an Eigen 3d Matrix for output
	SpiceDouble spice_C[3][3];
	pxform_c(from_frame.c_str(), to_frame.c_str(), et, spice_C);

	Eigen::Matrix3d C; //want to use Eigen vectors & matrices in the rest of the script
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			C(i, j) = spice_C[i][j];
		}
	}
	return C;
}
double SpiceHandler::str_date_to_et(std::string date_string)
{
	double et;
	str2et_c(date_string.c_str(), &et);

	return et;
}
#pragma endregion utilities


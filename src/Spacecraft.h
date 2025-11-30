
#pragma once
#include <string>
#include <vector>
#include <ostream>
#include "structure_definitions.h"
#include "Integrator.h"
#include "Planet.h"


class Spacecraft
{
public:
	Spacecraft(Integrator& integrator);
	Spacecraft(Integrator& integrator, std::string name, State state);
	Spacecraft(Integrator& integrator, std::string name, double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb);
	Spacecraft(Integrator& integrator, std::string name, double et, Eigen::Vector<double, 6> coes, double mu_cb);
	Spacecraft(const Spacecraft& other);
	Spacecraft& operator=(const Spacecraft& other);

	//getters
	std::string get_name() const;
	COE get_ref_conic() const;
	State get_state() const;
	const std::vector<double>& get_et_history() const;
	const std::vector<Eigen::Vector<double, 6>>& get_cartesian_history() const; //ICRF
	const std::vector<Eigen::Vector<double, 6>>& get_coe_history() const; //instantaneous elements

	//setters
	void set_name(std::string new_name);
	void set_state(State new_state); //sets the current_state and also updates the state_history vector

	void set_ref_conic(COE new_conic); //only to be used if the spacecraft is reset (may be reset to new orbit)
	void reset_state(State state0); //resets state_history to only the new state0 (and sets current_state accordingly)
	void reset_state(double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb); //alternative reset function for convenience
	void reset_state(double et, Eigen::Vector<double, 6> coes, double mu_cb);

	//utilities
	void update_current_state_coes(double mu_cb);
	void update_current_state_cart(double mu_cb);
	void reset_state_history_vecs(State new_state);
	void add_state_to_history_vecs(State new_state);
	//double elevation_to_ground_stations(Planet& planet); //inputs will be provided by constellation class
	void apply_dv(Eigen::Vector3d dv_vec);
	
	void step(double dt);

	void history_row_count_validation();
	void build_eigen_state_history();
	void write_history_to_csv(std::string filename);

private:

	std::string name;

	COE ref_conic; //orbital elements for reference, conic orbit assumed when initializing the satellite's state

	//note: for now, leaving mass out of this simulation as a first-pass, quick software example
	State current_state;
	std::vector<double> et_history;
	std::vector<Eigen::Vector<double, 6>> cartesian_history; //ICRF
	std::vector<Eigen::Vector<double, 6>> coe_history; //instantaneous orbital elements

	Eigen::MatrixXd collected_history;
	//note: will hold the state history in the et_history, cartesian_history, & coe_history vectors;
	//		this Eigen matrix will be built as needed for convenient vector math & data output

	Integrator& integrator;
};


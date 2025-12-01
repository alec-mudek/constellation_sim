
#pragma once
#include <string>
#include <vector>
#include <ostream>
#include <astrokit/math_utils.h>
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
	TrackingState get_tracking() const;
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

	std::size_t get_et_index(double et); //gets the location index of the closest match to a given et in the et_history vector
	void update_tracking(Spacecraft& neighbor1, Spacecraft& neighbor2); //fills in the tracking state information
	//note: only keeping current information for this atm, not storing the full history (by design)
	bool check_in_bounds(const BoundingBox& bounds); //returns true/false whether or not the spacecraft is currently within the bounds
	//note: returns an Eigen vector instead of a single bool so that we can know which check(s) failed

	//data handling
	void history_row_count_validation();
	Eigen::MatrixXd build_partial_eigen_history(std::size_t ix0, std::size_t ixf);
	void build_eigen_state_history();
	void write_history_to_csv(std::string filename);

private:

	std::string name;

	COE ref_conic; //orbital elements for reference, conic orbit assumed when initializing the satellite's state
	double ref_period; //period of the conic reference orbit

	//note: for now, leaving mass out of this simulation as a first-pass, quick software example
	State current_state;
	//note: using std::vector instead of Eigen matrices here since it's more efficient for memory reallocation as the states grow.
	//		will convert the states to Eigen matrices later as needed for efficient computations.
	std::vector<double> et_history;
	std::vector<Eigen::Vector<double, 6>> cartesian_history; //ICRF
	std::vector<Eigen::Vector<double, 6>> coe_history; //instantaneous orbital elements

	TrackingState tracking; //contains bounding box information for the current time

	Eigen::MatrixXd collected_history;
	//note: will hold the state history in the et_history, cartesian_history, & coe_history vectors;
	//		this Eigen matrix will be built as needed for convenient vector math & data output

	Integrator& integrator;
};


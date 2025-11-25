
#pragma once
#include <string>
#include <vector>
#include "structure_definitions.h"
#include "Integrator.h"


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
	State get_state() const;
	const std::vector<State>& get_state_history() const;

	//setters
	void set_name(std::string new_name);
	void set_state(State new_state); //sets the current_state and also updates the state_history vector

	void reset_state(State state0); //resets state_history to only the new state0 (and sets current_state accordingly)
	void reset_state(double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb); //alternative reset function for convenience
	void reset_state(double et, Eigen::Vector<double, 6> coes, double mu_cb);

	//utilities
	void update_current_state_coes(double mu_cb);
	void update_current_state_cart(double mu_cb);
	void apply_dv(Eigen::Vector3d dv_vec);
	void step(double dt);

private:

	std::string name;
	//note: for now, leaving mass out of this simulation as a first-pass, quick software example
	State current_state;
	std::vector<State> state_history;

	Integrator& integrator;
};


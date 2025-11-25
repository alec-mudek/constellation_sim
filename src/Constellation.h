
#pragma once
#include "Planet.h"
#include "Spacecraft.h"
#include "Integrator.h"

class Constellation
{
public:
	Constellation(Planet& cb, Integrator& integrator);
	Constellation(Planet& cb, Integrator& integrator, double et0);
	Constellation(Planet& cb, Integrator& integrator, double et0, std::vector<Spacecraft> sc_list);

	//going for a singleton-ish pattern for the Constellation class; don't want it to be copyable
	Constellation(const Constellation&) = delete;
	Constellation& operator=(const Constellation&) = delete;
	Constellation(Constellation&&) = delete;
	Constellation& operator=(Constellation&&) = delete;

	//getters
	const std::vector<Spacecraft>& get_sats() const;
	double get_et() const;

	//setters
	void set_et(double new_et);

	//utilities
	void add_spacecraft(Spacecraft new_sc);
	void add_spacecraft(std::string sc_name, State sc_state);
	void add_spacecraft(std::string sc_name, double et0, Eigen::Vector3d pos0, Eigen::Vector3d vel0);
	void add_spacecraft(std::string sc_name, double et0, Eigen::Vector<double, 6> coes);

	//normally I'd have a standalone propagation class but I'm just defining that code here
	// for the purposes of this project
	void propagate(double duration, double step_size);
	//note: want to propagate every spacecraft in the constellation for each step before moving on

private:
	Planet& cb;
	Integrator& integrator;

	std::vector<Spacecraft> spacecraft;
	double current_et;

};


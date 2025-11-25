#pragma once

#include "Planet.h"
#include "ForceModel.h"


class Integrator
{
public:
	Integrator(Planet& cb, ForceModel& fm);

	//going for a singleton-ish pattern for the Integrator class; don't want it to be copyable
	Integrator(const Integrator&) = delete;
	Integrator& operator=(const Integrator&) = delete;
	Integrator(Integrator&&) = delete;
	Integrator& operator=(Integrator&&) = delete;

	const Planet& get_cb() const;

	Eigen::Vector<double, 6> step(double t, double dt, const Eigen::Vector<double, 6>& state);

private:
	Planet& cb;
	ForceModel& fm;

};
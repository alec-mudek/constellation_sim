#pragma once

#include "Planet.h"	
#include <astrokit/force_models.h>
#include <astrokit/integrators.h>

class ForceModel
{
public:
	//normally I'd break integration up into it's own class but for simplicity's sake I'm just defining both the EOMs and 
	// integration step function here.
	//note: NOT set up for n-body perturbations as is; sticking with J2 only for this demo
	ForceModel(Planet& cb);
	ForceModel(Planet& cb, bool include_j2);

	//going for a singleton-ish pattern for the ForceModel class; don't want it to be copyable
	ForceModel(const ForceModel&) = delete;
	ForceModel& operator=(const ForceModel&) = delete;
	ForceModel(ForceModel&&) = delete;
	ForceModel& operator=(ForceModel&&) = delete;

	void set_include_j2(bool j2_included);

	Eigen::Vector<double, 6> eoms(double t, const Eigen::Vector<double, 6>& state);


private:
	Planet& cb; //central body;

	bool include_j2;

};

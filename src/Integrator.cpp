#include "Integrator.h"
#include <astrokit/integrators.h>


Integrator::Integrator(Planet& cb, ForceModel& fm) : cb(cb), fm(fm)
{
}

const Planet& Integrator::get_cb() const
{
	return this->cb;
}

Eigen::Vector<double, 6> Integrator::step(double t, double dt, const Eigen::Vector<double, 6>& state)
//fixed-step RK4
{
	//need to use a lambda to make ForceModel's EOMs method work with astrokit (needs to be callable)
	auto f = [this](double tt, const Eigen::Vector<double, 6>& yy)
	{
		return fm.eoms(tt, yy);
	};

	return astrokit::rk4_step(t, dt, state, f);
}

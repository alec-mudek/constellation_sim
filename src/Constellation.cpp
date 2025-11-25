#include "Constellation.h"
#include <astrokit/integrators.h>

Constellation::Constellation(Planet& cb, Integrator& integrator) : cb(cb), current_et(0.0), spacecraft{}, integrator(integrator)
{
}

Constellation::Constellation(Planet& cb, Integrator& integrator, double et0) : cb(cb), integrator(integrator)
{
	set_et(et0);
}

Constellation::Constellation(Planet& cb, Integrator& integrator, double et0, std::vector<Spacecraft> sc_list) : cb(cb), integrator(integrator)
{
	set_et(et0);
	this->spacecraft = sc_list; //provided a vector of already-initialized s/c to the constructor
}

#pragma region getters
const std::vector<Spacecraft>& Constellation::get_sats() const
{
	return this->spacecraft;
}

double Constellation::get_et() const
{
	return this->current_et;
}
#pragma endregion getters

#pragma region setters
void Constellation::set_et(double new_et)
{
	this->current_et = new_et;
}
#pragma endregion setters

#pragma region utilities
void Constellation::add_spacecraft(Spacecraft new_sc)
{
	this->spacecraft.push_back(new_sc);
}

void Constellation::add_spacecraft(std::string sc_name, State sc_state)
{
	this->spacecraft.push_back(Spacecraft(this->integrator, sc_name, sc_state));
}

void Constellation::add_spacecraft(std::string sc_name, double et0, Eigen::Vector3d pos0, Eigen::Vector3d vel0)
{
	this->spacecraft.push_back(Spacecraft(this->integrator, sc_name, et0, pos0, vel0, cb.get_mu()));
}

void Constellation::add_spacecraft(std::string sc_name, double et0, Eigen::Vector<double, 6> coes)
{
	this->spacecraft.push_back(Spacecraft(this->integrator, sc_name, et0, coes, cb.get_mu()));
}

void Constellation::propagate(double duration, double step_size)
{
	//the Spacecraft class has its own Step() function which is responsible for applying an rk4 step to the 
	// current state and updating the Spacecraft's state_history accordingly
	//the Constellation class is responsible for looping through it's Spacecraft vector and calling each 
	// Spacecraft's Step() method for each time step.
	//note: it would be easier to just fully propagate each spacecraft at a time then compare time histories,
	//      but propagating the full constellation together gives more modeling flexibility for future features
	//also note: for now, we're assuming forward propagation only here
	double total_time = 0.0;
	while (total_time + step_size < duration)
	{
		for (auto& sc : this->spacecraft)
		{
			sc.step(step_size);
		}
		total_time += step_size;
	}
	if (total_time < duration) //we need one more partial step
	{
		double partial_step = duration - total_time;
		for (auto& sc : this->spacecraft)
		{
			sc.step(partial_step);
		}
	}
}
#pragma endregion utilities



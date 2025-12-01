#include "Constellation.h"
#include <astrokit/integrators.h>

Constellation::Constellation(Planet& cb, Integrator& integrator) : 
	cb(cb), current_et(0.0), spacecraft{}, sc_bounds{}, integrator(integrator)
{
}

Constellation::Constellation(Planet& cb, Integrator& integrator, double et0) : 
	cb(cb), spacecraft{}, sc_bounds{}, integrator(integrator)
{
	set_et(et0);
}

Constellation::Constellation(Planet& cb, Integrator& integrator, double et0, std::vector<Spacecraft> sc_list, BoundingBox sc_bounds) : 
	cb(cb), integrator(integrator)
{
	set_et(et0);
	this->spacecraft = sc_list; //provided a vector of already-initialized s/c to the constructor
	set_sc_bounds(sc_bounds); //provided the mean orbital element bounds for each spacecraft
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

BoundingBox Constellation::get_sc_bounds() const
{
	return this->sc_bounds;
}
#pragma endregion getters

#pragma region setters
void Constellation::set_et(double new_et)
{
	this->current_et = new_et;
}

void Constellation::set_sc_bounds(BoundingBox new_bounds)
{
	this->sc_bounds = new_bounds;
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
	if (total_time < duration) //we need one more partial step; want to always include the exact final time in the output
	{
		double partial_step = duration - total_time;
		for (auto& sc : this->spacecraft)
		{
			sc.step(partial_step);
		}
	}
}
void Constellation::save_spacecraft_histories(std::string file_name_root)
{
	//store spacecraft names as they're used to name the output files; want to make sure there are no duplicates
	std::vector<std::string> used_names{};

	//loop through each spacecraft in the constellation and have them write their state histories
	for (auto& sc : this->spacecraft)
	{
		std::string sc_name = sc.get_name();

		//check how many times the spacecraft name has already appeared
		int counter = 0;
		for (std::size_t i = 0; i < used_names.size(); i++)
		{
			if (sc_name == used_names[i])
			{
				counter++;
			}
		}
		std::string file_name = file_name_root + sc_name;
		if (counter > 0)
		{
			file_name += std::to_string(counter);
		}
		sc.write_history_to_csv(file_name + ".csv"); //appending the csv file type here instead of in Spacecraft; may change later
	}
}
#pragma endregion utilities



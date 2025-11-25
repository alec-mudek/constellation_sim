#include "Spacecraft.h"
#include <astrokit/state_converter.h>

Spacecraft::Spacecraft(Integrator& integrator) : name("Default"), current_state{}, state_history(), integrator(integrator)
{	
}

//best option is to provide the State struct directly in the constructor
Spacecraft::Spacecraft(Integrator& integrator, std::string name, State state) : integrator(integrator)
{
	set_name(name);
	reset_state(state); //ResetState function sets the current_state and stores it as the first (and only) entry in the state_history
}

//there may be times it is convenient to just provide the cartesian state (& mu) and let the constructor fill in the COE information
Spacecraft::Spacecraft(Integrator& integrator, std::string name, double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb) :
	integrator(integrator)
{
	set_name(name);
	reset_state(et, pos, vel, mu_cb);
}

//there will also be times we want to initialize a spacecraft by COEs
Spacecraft::Spacecraft(Integrator& integrator, std::string name, double et, Eigen::Vector<double, 6> coes, double mu_cb) :
	integrator(integrator)
{
	set_name(name);
	reset_state(et, coes, mu_cb);
}

Spacecraft::Spacecraft(const Spacecraft& other)
	: name(other.name),
	current_state(other.current_state),
	state_history(other.state_history),
	integrator(other.integrator)   // bind our reference to the same Integrator; note: the integrator handling is what requires these
{
}

Spacecraft& Spacecraft::operator=(const Spacecraft& other)
{
	if (this != &other)
	{
		this->name = other.name;
		this->current_state = other.current_state;
		this->state_history = other.state_history;
	}
	return *this;
}


#pragma region getters
std::string Spacecraft::get_name() const
{
	return this->name;
}

State Spacecraft::get_state() const
{
	return this->current_state;
}

const std::vector<State>& Spacecraft::get_state_history() const
{
	return this->state_history;
}
#pragma endregion getters

#pragma region setters
void Spacecraft::set_name(std::string new_name)
{
	this->name = new_name;
}

void Spacecraft::set_state(State new_state)
{
	//method for updating the spacecraft's state information
	//predominantly used by the propagator to update the s/c with each step
	this->current_state = new_state;
	this->state_history.push_back(new_state);
}

void Spacecraft::reset_state(State state0)
{
	this->current_state = state0;
	this->state_history = { state0 };
}

void Spacecraft::reset_state(double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb)
{
	//alternative ResetState method that doesn't require COEs to be provided
	this->current_state = State{ et, pos, vel, 0, 0, 0, 0, 0, 0 }; //will fill the COEs in with the helper function UpdateCurrentStateCOEs
	update_current_state_coes(mu_cb); //fills the COE information in in place

	this->state_history = { this->current_state };
	
}

void Spacecraft::reset_state(double et, Eigen::Vector<double, 6> coes, double mu_cb)
{
	double sma = coes[0];
	double ecc = coes[1];
	double inc = coes[2];
	double raan = coes[3];
	double argp = coes[4];
	double ta = coes[5];

	this->current_state = State{ et, Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, sma, ecc, inc, raan, argp, ta };
	update_current_state_cart(mu_cb);

	this->state_history = { this->current_state };
}
#pragma endregion setters

#pragma region utilities
void Spacecraft::update_current_state_coes(double mu_cb)
{
	//function to take the current_state and recalculate the COEs based on the position and velocity vectors
	// note: often make changes to the cartesian state, then compute the COEs from the new pos, vel

	//use astrokit's Cart2COE function
	Eigen::Vector<double, 6> cart;
	cart << this->current_state.pos, this->current_state.vel;

	Eigen::Vector<double, 6> coes = astrokit::cart_to_coe(cart, mu_cb);
	this->current_state.sma = coes[0];
	this->current_state.ecc = coes[1];
	this->current_state.inc = coes[2];
	this->current_state.raan = coes[3];
	this->current_state.argp = coes[4];
	this->current_state.ta = coes[5];
}

void Spacecraft::update_current_state_cart(double mu_cb)
{
	Eigen::Vector<double, 6> coes;
	coes << this->current_state.sma, this->current_state.ecc, this->current_state.inc,
			this->current_state.raan, this->current_state.argp, this->current_state.ta;

	Eigen::Vector<double, 6> cart = astrokit::coe_to_cart(coes, mu_cb);
	this->current_state.pos = cart.segment<3>(0);
	this->current_state.vel = cart.segment<3>(3);
}

void Spacecraft::apply_dv(Eigen::Vector3d dv_vec)
{
	// note: applies an impulsive delta-v to the sc
	// need to update both the velocity vector and COEs
	//don't want to require mu as an input, so reverse engineer it before updating the velocity
	double R = this->current_state.pos.norm();
	double V = this->current_state.vel.norm();
	double mu_cb = V * V / (2 / R - 1 / this->current_state.sma);

	//now update the velocity vector
	this->current_state.vel = this->current_state.vel + dv_vec;

	//and update the coes
	update_current_state_coes(mu_cb);

	//now that our state is fully updated; append it as a new step in our state_history
	this->state_history.push_back(this->current_state);
}
void Spacecraft::step(double dt)
{
	double t = this->current_state.et;
	Eigen::Vector<double, 6> state;
	state << this->current_state.pos, this->current_state.vel;

	Eigen::Vector<double, 6> new_state = integrator.step(t, dt, state);

	//update the time & cartesian components of the current_state
	this->current_state.et += dt;
	this->current_state.pos = new_state.segment<3>(0);
	this->current_state.vel = new_state.segment<3>(3);
	//and fill in the new COE information
	update_current_state_coes(this->integrator.get_cb().get_mu()); 

	//now add the updated state to the state_history
	this->state_history.push_back(this->current_state);
}
#pragma endregion utilities
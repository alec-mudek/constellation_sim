#include <fstream>
#include <stdexcept>
#include "Spacecraft.h"
#include <astrokit/state_converter.h>

Spacecraft::Spacecraft(Integrator& integrator) : 
	name("Default"), current_state{}, et_history(), cartesian_history(), coe_history(), integrator(integrator)
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
	et_history(other.et_history),
	cartesian_history(other.cartesian_history),
	coe_history(other.coe_history),
	integrator(other.integrator)   // bind our reference to the same Integrator; note: the integrator handling is what requires these
{
}

Spacecraft& Spacecraft::operator=(const Spacecraft& other)
{
	if (this != &other)
	{
		this->name = other.name;
		this->current_state = other.current_state;
		this->et_history = other.et_history;
		this->cartesian_history = other.cartesian_history;
		this->coe_history = other.coe_history;
	}
	return *this;
}


#pragma region getters
std::string Spacecraft::get_name() const
{
	return this->name;
}

COE Spacecraft::get_ref_conic() const
{
	return this->ref_conic;
}

State Spacecraft::get_state() const
{
	return this->current_state;
}

const std::vector<double>& Spacecraft::get_et_history() const
{
	return this->et_history;
}

const std::vector<Eigen::Vector<double, 6>>& Spacecraft::get_cartesian_history() const
{
	return this->cartesian_history;
}

const std::vector<Eigen::Vector<double, 6>>& Spacecraft::get_coe_history() const
{
	return this->coe_history;
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
	add_state_to_history_vecs(new_state);
}

void Spacecraft::set_ref_conic(COE new_conic)
{
	this->ref_conic = new_conic;
}

void Spacecraft::reset_state(State state0)
{
	//set the current state
	this->current_state = state0;

	//and reinitialize the state history vectors
	reset_state_history_vecs(state0);

	//and set the reference conic
	COE reference(state0.sma, state0.ecc, state0.inc, state0.raan, state0.argp, state0.ta);
	set_ref_conic(reference);
}

void Spacecraft::reset_state(double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb)
{
	//alternative ResetState method that doesn't require COEs to be provided
	this->current_state = State{ et, pos, vel, 0, 0, 0, 0, 0, 0 }; //will fill the COEs in with the helper function UpdateCurrentStateCOEs
	update_current_state_coes(mu_cb); //fills the COE information in in place

	reset_state_history_vecs(this->current_state);

	COE reference(get_state().sma, get_state().ecc, get_state().inc, get_state().raan, get_state().argp, get_state().ta);
	set_ref_conic(reference);
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

	reset_state_history_vecs(this->current_state);

	COE reference(get_state().sma, get_state().ecc, get_state().inc, get_state().raan, get_state().argp, get_state().ta);
	set_ref_conic(reference);
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

void Spacecraft::reset_state_history_vecs(State new_state)
{
	//et history for time
	this->et_history = { new_state.et };

	//cartesian for ICRF pos, vel
	Eigen::Vector<double, 6> cart_state;
	cart_state << new_state.pos, new_state.vel;
	this->cartesian_history = { cart_state };

	//coe for instaneous orbital elements
	Eigen::Vector<double, 6> coe_state;
	coe_state << new_state.sma, new_state.ecc, new_state.inc, new_state.raan, new_state.argp, new_state.ta;
	this->coe_history = { coe_state };
}

void Spacecraft::add_state_to_history_vecs(State new_state)
{
	this->et_history.push_back(new_state.et);

	Eigen::Vector<double, 6> cart_state;
	cart_state << new_state.pos, new_state.vel;
	this->cartesian_history.push_back(cart_state);

	Eigen::Vector<double, 6> coe_state;
	coe_state << new_state.sma, new_state.ecc, new_state.inc, new_state.raan, new_state.argp, new_state.ta;
	this->coe_history.push_back(coe_state);
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
	add_state_to_history_vecs(this->current_state);
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
	add_state_to_history_vecs(this->current_state);
}

void Spacecraft::history_row_count_validation()
{
	//want to throw an error if the state histories somehow ended up as different sizes
	if (this->et_history.size() != this->cartesian_history.size() || this->et_history.size() != this->coe_history.size())
	{
		throw std::runtime_error("History length mismatch in Spacecraft " + get_name() + ".");
	} //shouldn't be necessary but a good sanity check during development
}

void Spacecraft::build_eigen_state_history()
{
	history_row_count_validation(); //sanity check

	//before we store the data, need to resize the collected_history matrix
	std::size_t n_states = get_et_history().size();
	this->collected_history.resize(n_states, 13);

	//now store the vector data into the Eigen matrix
	for (std::size_t i = 0; i < n_states; i++)
	{
		this->collected_history(i, 0) = this->et_history[i];
		this->collected_history.block<1, 6>(i, 1) = this->cartesian_history[i].transpose();
		this->collected_history.block<1, 6>(i, 7) = this->coe_history[i].transpose();
	}
}

void Spacecraft::write_history_to_csv(std::string filename)
{
	history_row_count_validation(); //sanity check

	//make sure the collected_history is up-to-date
	build_eigen_state_history();
	
	//now create a csv file to save to
	std::ofstream f(filename);

	//write the header row
	f << "et,rx,ry,rz,vx,vy,vz,sma,ecc,inc,raan,argp,ta\n";

	//and write the data
	Eigen::IOFormat csv(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");
	f << this->collected_history.format(csv);
}
#pragma endregion utilities
#include <fstream>
#include <stdexcept>
#include "Spacecraft.h"
#include <astrokit/state_converter.h>

Spacecraft::Spacecraft(Integrator& integrator) : 
	name("Default"), current_state{}, et_history(), cartesian_history(), coe_history(), integrator(integrator)
{	
}

//best option is to provide the State struct directly in the constructor
Spacecraft::Spacecraft(Integrator& integrator, std::string name, State state) : 
	integrator(integrator)
{
	set_name(name);
	reset_state(state); //reset_state function sets the current_state and stores it as the first (and only) entry in the state_history
	//reset_state also handles reseting the current tracking state
}

//there may be times it is convenient to just provide the cartesian state (& mu) and let the constructor fill in the COE information
Spacecraft::Spacecraft(Integrator& integrator, std::string name, double et, Eigen::Vector3d pos, Eigen::Vector3d vel, double mu_cb) :
	integrator(integrator)
{
	set_name(name);
	reset_state(et, pos, vel, mu_cb); //overloaded functions handle necessary computations to fill in the rest of the State
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
	tracking(other.tracking),
	collected_history(other.collected_history),
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
		this->tracking = other.tracking;
		this->collected_history = other.collected_history;
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

TrackingState Spacecraft::get_tracking() const
{
	return this->tracking;
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
	this->ref_period = 2 * astrokit::PI * sqrt(pow(new_conic.sma, 3)  / this->integrator.get_cb().get_mu());
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

std::size_t Spacecraft::get_et_index(double target_et)
{
	//the main use of this function is to find the vector indeces for orbital element averages.
	//  to that end, we don't need an exact match. just the closest points in our dataset. if the
	//  exact target_et doesn't appear in the et_history, this function will loop until 
	
	//before we loop, make sure the target_et is in the et_history range (ASSUMES CONSISTENT FORWARD/BACKPROP DIRECTION FOR HISTORY)
	double et0 = this->et_history[0];
	double etf = this->et_history[this->et_history.size() - 1];

	//check that they haven't asked for a future time step in the propagation (either forward in time for forward prop or backward for backprop)
	if ((etf > et0 && target_et > etf) || (et0 > etf && target_et < etf))
	{
		throw std::runtime_error("Spacecraft " + get_name() + " requested state beyond propagation bounds.");
	}

	//on the other hand, if they requested a target_et that is in the right direction but extends beyond the state history
	if ((etf > et0 && target_et < et0) || (et0 > etf && target_et > et0))
	{
		return 0; //just return index 0 (just use the full history in this case)
	}

	//if we make it here, the target_et should fall somewhere in our propagation history
	//loop from the end to the beginning of et_history and find the closest match
	int current_sign = astrokit::sign(target_et - this->current_state.et);
	int previous_sign;
	for (std::size_t i = this->et_history.size()-1; i >= 0; i--)
	{
		//if not, check for if we've gone from undershooting to overshooting our target_et
		previous_sign = current_sign;
		current_sign = astrokit::sign(target_et - this->et_history[i]);
		//note: the astrokit::sign function returns +1, 0, or -1 if the input is >0, =0, or <0 respectively
		if (current_sign != previous_sign)
		{
			return i;
		}
	}
	//should never make it here, but just in case (REMOVE BEFORE FLIGHT)
	throw std::runtime_error("Error in Spacecraft " + get_name() + " get_et_index function for time " + std::to_string(target_et) + ".");
}

void Spacecraft::update_tracking(Spacecraft& neighbor1, Spacecraft& neighbor2)
{
	//first need to determine the mean elements over the last orbit period
	double etf = get_state().et;
	double et0 = etf - this->ref_period;

	std::size_t start_ix = get_et_index(et0);
	std::size_t stop_ix = get_et_index(etf);

	//get the data we need for the average elements
	//note: computing the mean over the last period's worth of data (using the period of the reference conic)
	Eigen::MatrixXd dat = build_partial_eigen_history(start_ix, stop_ix);
	
	//now fill in the tracking data
	//dat columns:
	//et, rx, ry, rz, vx, vy, vz, sma, ecc, inc, raan, argp, ta
	// 0,  1,  2,  3,  4,  5,  6,   7,   8,   9,   10,   11, 12
	this->tracking.et = etf;
	this->tracking.sma_mean = dat.col(7).mean();
	this->tracking.inc_mean = dat.col(9).mean();
	this->tracking.raan_mean = dat.col(10).mean();

	//now check phasing w/neighbors
	//note: don't want to worry about angle wrapping issues; just use the position vectors and find the angle between
	Eigen::MatrixXd n1_dat = neighbor1.build_partial_eigen_history(start_ix, stop_ix);
	Eigen::MatrixXd n2_dat = neighbor2.build_partial_eigen_history(start_ix, stop_ix);
	double n1_angle_sum = 0.0;
	double n2_angle_sum = 0.0;
	for (std::size_t i = 0; i < dat.size(); i++)
	{
		Eigen::Vector3d my_pos = dat.row(i).segment<3>(1); //segment of 3 columns starting at index 1 in row i -> position components at the correct time step
		double n1_angle = astrokit::angle_between_vecs(my_pos, n1_dat.row(i).segment<3>(1));
		double n2_angle = astrokit::angle_between_vecs(my_pos, n2_dat.row(i).segment<3>(1));

		n1_angle_sum += n1_angle;
		n2_angle_sum += n2_angle;
	}
	this->tracking.neighbor1_rel_angle = n1_angle_sum / dat.size();
	this->tracking.neighbor2_rel_angle = n2_angle_sum / dat.size();
}

bool Spacecraft::check_in_bounds(const BoundingBox& bounds)
{
	//working on this logic now
	//using Galileo as a reference point: GALILEO CONSTELLATION: EVALUATION OF STATION KEEPING STRATEGIES by Navarro-Reyes, et. al.
	return true;
}
#pragma endregion utilities

#pragma region data handling
void Spacecraft::history_row_count_validation()
{
	//want to throw an error if the state histories somehow ended up as different sizes
	if (this->et_history.size() != this->cartesian_history.size() || this->et_history.size() != this->coe_history.size())
	{
		throw std::runtime_error("History length mismatch in Spacecraft " + get_name() + ".");
	} //shouldn't be necessary but a good sanity check during development
}

Eigen::MatrixXd Spacecraft::build_partial_eigen_history(std::size_t ix0, std::size_t ixf)
{
	const int n_rows = ixf - ix0;
	Eigen::MatrixXd out(ixf - ix0, 13);

	for (std::size_t i = ix0; i <= ixf; i++)
	{
		out(i, 0) = this->et_history[i];
		out.block<1, 6>(i, 1) = this->cartesian_history[i].transpose();
		out.block<1, 6>(i, 7) = this->coe_history[i].transpose();
	}
	return out;
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
#pragma endregion data handling
/*
Script to intialize and propagate a Walker Delta Earth Constellation

Current modeling assumptions:
-Perfect spheroid Earth for ground station surface normal, ground tracks, etc.
-Two-body propagation + J2 perturbations for the force model
-Purely an elevation contraint for ground station access

For now, simulation includes a single ground station and a 12/3/1 Walker Delta constellation.
Picked 56 deg inclination and 25,000 km sma.
*/

#include "Planet.h"
#include "SpiceHandler.h"
#include "GroundStation.h"
#include "WalkerDelta.h"
#include <astrokit/constants.h>

#include <iostream>

int main()
{
	SpiceHandler spice;

	//time info
	std::string start_date = "Jan 01 2026 00:00:000";
	double et0 = spice.str_date_to_et(start_date); //convert start_date string to spice ephemeris time
	double prop_time = 86400*3; //[s]; 3 day prop
	double prop_step = 10; //[s]

	//object initialization
	Planet earth(spice, astrokit::EARTH.MU_km3_s2, astrokit::EARTH.R_MEAN_km, astrokit::EARTH.R_EQUATOR_km, astrokit::EARTH.J2, 399, "IAU_EARTH");
	ForceModel fm(earth);
	Integrator rk4(earth, fm);
	WalkerDelta wd_const(earth, rk4, et0, 27, 3, 1, 56.0 * astrokit::DEG2RAD, 29600.0); //inc & sma roughly based off Galileo WD constellation

	//define ground station(s)
	//for now, just using approxiamte APL lon/lat
	earth.new_station("APL", 39.0 * astrokit::DEG2RAD, -77.0 * astrokit::DEG2RAD);

	//now have both our ground station and constellation satellites initialized
	//time to propagate
	wd_const.propagate(prop_time, prop_step);

	std::cout << "Propagation finished.";
	//NEXT STEP: Data processing, analysis, & visualization

	//output csv state histories for each satellite
	wd_const.save_spacecraft_histories("C:/constellation_sim_results/");

	return 0;
}
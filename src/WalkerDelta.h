#pragma once

#include "Constellation.h"

class WalkerDelta : public Constellation
{
public:
	WalkerDelta(Planet& cb, Integrator& integrator, double et0, int T, int P, int F, double inc, double sma);
	WalkerDelta(Planet& cb, Integrator& integrator, double et0, int T, int P, int F, double inc, double sma, double raan0); //allows specifying the reference raan (defaults to 0)
	//a specific type of constellation which dictates how the satellites are initialized
	//only really need a constructor for this, the constellation logic itself is all the same as the Constellation class
	//breaking the satellite intialization out into its own method though

	void initialize_satellites(int number_of_sats, int number_of_planes, int phasing, double inc, double sma, double raan0);
	// Walker Delta T/P/F Convention
	//T = number of satellites (dispersed evenly across true anomly & orbital planes)
	//P = number of orbital planes (dispersed evenly across raan)
	//F = phasing parameter b/w orbital planes (true anomaly offset)
};
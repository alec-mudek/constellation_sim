#include "walkerdelta.h"
#include "Spacecraft.h"
#include <astrokit/constants.h>

WalkerDelta::WalkerDelta(Planet& cb, Integrator& integrator, double et0, int T, int P, int F, double inc, double sma) :
	Constellation(cb, integrator, et0)
{
	initialize_satellites(T, P, F, inc, sma, 0.0); //default to raan0 = 0 if unspecified
}

WalkerDelta::WalkerDelta(Planet& cb, Integrator& integrator, double et0, int T, int P, int F, double inc, double sma, double raan0) :
	Constellation(cb, integrator, et0)
{
	initialize_satellites(T, P, F, inc, sma, raan0);
}

void WalkerDelta::initialize_satellites(int number_of_sats, int number_of_planes, int phasing, double inc, double sma, double raan0)
{
	//Walker Delta constellation consists of satellites in circular orbits evenly distributed in RAAN & 
	// true anomaly (phased across orbital planes). All planes share the same inclination & sma.
	// 
	//T = number of satellites (dispersed evenly across true anomly & orbital planes)
	//P = number of orbital planes (dispersed evenly across raan)
	//F = phasing parameter b/w orbital planes (true anomaly offset)
	//inc = inclination for all orbits
	//sma = semi-major axis for all orbits
	//raan0 = reference RAAN; used for the first orbit plane

	assert(number_of_sats % number_of_planes == 0 && "WalkerDelta requires T (number_of_sats) to be divisible by P (number_of_planes)");
	//note: not building in "what to do with leftover spacecraft?" logic for now

	int sats_per_plane = number_of_sats / number_of_planes;
	double raan_spacing = 2.0 * astrokit::PI / number_of_planes;
	double ta_spacing = 2.0 * astrokit::PI / sats_per_plane;
	//note: sticking with "true anomaly" for consistency with generalized COE convention; this is
	//		really the argument of latitude since we're in circular orbits (measured from RAAN crossing)
	double ta0_offset = 2.0 * astrokit::PI * phasing / number_of_sats;
	//note: the phasing term provides a fractional offset for the reference true anomaly (the ta of the
	//		first spacecraft) in each plane. 
	//      For example, if you create a walker delta constellation with 12 satellites, then an F value 
	//		of 1 would offset each reference true anomly by 1/12 of an orbit (2pi*1/12). The first 
	//		spacecraft in the first orbital plane may start at ta=0. The first spacecraft in the second
	//		orbital plane would start at ta=pi/6. The first spacecraft in the third orbital plane would
	//		start at ta=2pi/6. Etc...
	
	//now create each satellite; will use COEs to define their initial state
	// sma & inc are fixed; by definition argp=e=0
	// need to find raan & ta for each satellite
	//loop through the orbital planes
	for (int plane_number = 0; plane_number < number_of_planes; plane_number++)
	{
		double raan = raan0 + raan_spacing * plane_number; //all satellites in each plane will share the same raan
		double ta0 = ta0_offset * plane_number;
		for (int sat_number = 0; sat_number < sats_per_plane; sat_number++)
		{
			std::string sat_name = "plane_" + std::to_string(plane_number) + "_sat_" + std::to_string(sat_number);
			double ta = ta0 + ta_spacing * sat_number; //true anomaly/argument of latitude is what provides unique state to each sat
			
			Eigen::Vector<double, 6> coes;
			coes << sma, 0, inc, raan, 0, ta;

			add_spacecraft(sat_name, get_et(), coes);
		}
	}
}

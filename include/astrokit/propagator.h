#pragma once

#include <vector>
#include <cmath>

namespace astrokit
{
	template <typename State>
	struct StateHistory
	{
		std::vector<double> t;
		std::vector<State> y;
	};

	template <typename State, typename StepMethod>
	inline StateHistory<State> propagator(double t0, double dt, double tf, const State& y0, const StepMethod& step)
	{
		///intialize state_history
		StateHistory<State> state_history;

		//start with initial state
		double t = t0;
		State y = y0;

		state_history.t.push_back(t0);
		state_history.y.push_back(y0);
		
		//clarify prop direction; need to worry about forward vs. backward prop
		//e.g. if user says t0 = 100, tf = 0, dt = 5 we need to make sure to correctly parse that to steps of -5)
		const double prop_dir = (t0 <= tf) ? 1.0 : -1.0;
		dt = prop_dir * std::abs(dt);

		//now ready to start stepping
		while ((prop_dir > 0.0 && t + dt < tf) || (prop_dir < 0.0 && t + dt > tf))
		{
			y = step(t, dt, y);
			t += dt;

			state_history.t.push_back(t);
			state_history.y.push_back(y);
		}
		
		//unless the time span is perfectly divisible by dt, we'll have a remainder
		double dt_final = tf - t;

		if (dt_final != 0.0)
		{
			y = step(t, dt_final, y);
			t += dt_final;

			state_history.t.push_back(t);
			state_history.y.push_back(y);
		}

		return state_history;
	}
} // namespace astrokit
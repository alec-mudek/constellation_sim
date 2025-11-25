#pragma once

namespace astrokit 
{
	template <typename State, typename F>
	inline State rk4_step(double t, double dt, const State& y, const F& f)
	//note: just following the standard butcher tableau here
	//also note: this template will accept any State type but it must allow for addition and scalar multiplication
	//			 i.e. you can't use a std::array with this but you can use Eigen::Vector or EMTG::math::Matrix
	//also also note: anywhere this rk4_step function is used, the integration function, f, must be of the form f(t, y)
	{
		State k1 = f(t, y);
		State k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1);
		State k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2);
		State k4 = f(t + dt, y + dt * k3);

		return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
} // namespace astrokit

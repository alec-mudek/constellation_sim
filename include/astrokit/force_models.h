#pragma once

#include <Eigen/Dense>
#include <cmath>

namespace astrokit 
{
    //note: The following accel functions can't be used directly with the rk4_step function in integrators.h
    //      since they aren't of the form f(t, y). The idea is these will be wrapped in an integration method 
    //      in whatever external script they are used (more modular & generic this way).

    //basic, two-body motion
    inline Eigen::Vector<double, 6> accel_kep(const Eigen::Vector<double, 6>& state, double mu)
    {
        Eigen::Vector3d r = state.segment<3>(0);
        Eigen::Vector3d v = state.segment<3>(3);

        const double R = r.norm();
        const double R3 = R * R * R;

        Eigen::Vector<double, 6> dstate_dt;
        dstate_dt << v, -mu * r / R3 ;

        return dstate_dt;
    }

    // J2 gravity perturbation for oblate spheroid, z-axis aligned with spin axis
    // note: ONLY includes the J2 perturbation; no central body gravity included
    inline Eigen::Vector<double, 6> accel_j2(const Eigen::Vector<double, 6>& state, double mu, double Re, double J2)
    {
        Eigen::Vector3d r = state.segment<3>(0);
        Eigen::Vector3d v = state.segment<3>(3);

        const double R1 = r.norm();
        const double R2 = R1 * R1;
        const double R5 = R2 * R2 * R1;
        const double z2 = r(2) * r(2);

        const double factor = 1.5 * J2 * mu * Re * Re / R5;
        const double k = 5.0 * z2 / R2;

        Eigen::Vector3d a((r(0)) * (k - 1.0) * factor, (r(1)) * (k - 1.0) * factor, (r(2)) * (k - 3.0) * factor);

        Eigen::Vector<double, 6> dstate_dt;
        dstate_dt <<  v, a ;

        return dstate_dt;
    }

} // namespace astrokit

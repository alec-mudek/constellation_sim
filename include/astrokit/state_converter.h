#pragma once

#include <cmath>
#include <Eigen/Dense>
#include "math_utils.h"
#include "rotations.h"

namespace astrokit
{
    inline Eigen::Vector<double, 6> cart_to_coe(const Eigen::Vector<double, 6>& cart, const double mu)
    // slightly modified version of EMTG's StateConverter CartesiantoCOE function
    // return order: sma, ecc, inc, raan, argp, ta
    {
        Eigen::Vector3d r = cart.segment<3>(0);
        Eigen::Vector3d v = cart.segment<3>(3);

        double R = r.norm();
        double V = v.norm();

        // sma - use the orbital specific energy
        double eps = V * V / 2 - mu / R;
        double sma = -mu / (2 * eps);

        // ecc - directly compute ecc vector
        Eigen::Vector3d h = r.cross(v);
        Eigen::Vector3d e_vec = v.cross(h) / mu - r / R;
        double ecc = e_vec.norm();

        // inc - from angular momentum vector
        Eigen::Vector3d h_hat = h.normalized();
        double inc = safe_acos(h_hat[2]);

        // raan - compute ascending node vector
        // argp - use node & ecc vectors
        // ta - use position & eccentricity vector or n vector & eccentricity vector
        Eigen::Vector3d n = k_hat.cross(h_hat); //k_hat from rotations.h
        double N = n.norm(); //if equatorial orbit, k_hat = h_hat, raan = 0, and the math can breakdown
        double raan = 0; //default to zero; fill in if able
        double argp = 0;
        double ta   = 0;
        double delta = 1e-6;
        //is it non-equatorial?
        if (N > delta)
        {
            raan = safe_acos(n[0] / N);
            
            if (n[1] < 0.0)
            {
                raan = 2.0 * PI - raan;
            }

            //is it non-circular?
            if (ecc > delta)
            {
                argp = safe_acos(n.dot(e_vec) / (N * ecc));
                if (e_vec[2] < 0.0)
                {
                    argp = 2.0 * PI - argp;
                }
                ta = safe_acos(e_vec.dot(r) / (ecc * R));
                if (r.dot(v) < 0.0) //quadrant check
                {
                    ta = 2.0 * PI - ta;
                }
            }
            else //it's circular -> use the ascending node to define ta
            {
                ta = safe_acos(n.dot(r) / (N * R));
                if (r[2] < 0.0) //quadrant check
                {
                    ta = 2.0 * PI - ta;
                }
            }

        }
        else //trajectory is equatorial
        {
            //is it non-circular
            if (ecc > delta) //equatorial & non-circ
            {
                argp = safe_acos(e_vec[0] / ecc);
                if (e_vec[1] < 0.0) //quad check
                {
                    argp = 2.0 * PI - argp;
                }
                ta = safe_acos(e_vec.dot(r) / (ecc * R));
                if (r.dot(v) < 0.0) //quad check
                {
                    ta = 2.0 * PI - ta;
                }
            }
            else //equatorial, circular orbit
            {
                ta = safe_acos(r[2] / R);
                if (r[2] < 0.0) //quad check
                {
                    ta = 2.0 * PI - ta;
                }
            }
        }
        //wrap ta if needed
        if (std::abs(ta - 2.0 * PI) < delta)
        {
            ta = 0.0;
        }

        Eigen::Vector<double, 6> coes;
        coes << sma, ecc, inc, raan, argp, ta;

        return coes;
    }

    inline Eigen::Vector<double, 6> coe_to_cart(const Eigen::Vector<double, 6>& coes, const double mu)
    {
        double a = coes[0];
        double e = coes[1];
        double i = coes[2];
        double Om = coes[3];
        double w = coes[4];
        double f = coes[5];

        //start with perifocal coords
        double p = a * (1 - e * e);
        double r = p / (1 + e * std::cos(f));
        Eigen::Vector3d r_pf(r * std::cos(f), r * std::sin(f), 0);
        Eigen::Vector3d v_pf(-std::sin(f), e + std::cos(f), 0);
        v_pf *= std::sqrt(mu / p);

        //and construct the 3-1-3 rotation matrix
        double cOm = std::cos(Om);
        double sOm = std::sin(Om);
        double cw = std::cos(w);
        double sw = std::sin(w);
        double ci = std::cos(i);
        double si = std::sin(i);

        Eigen::Matrix3d C;
        C << cOm * cw - sOm * sw * ci, -cOm * sw - sOm * cw * ci, sOm* si,
             sOm* cw + cOm * sw * ci, -sOm * sw + cOm * cw * ci, -cOm * si,
             sw* si, cw* si, ci;

        Eigen::Vector3d r_cart = C * r_pf;
        Eigen::Vector3d v_cart = C * v_pf;

        Eigen::Vector<double, 6> cart;
        cart << r_cart, v_cart;

        return cart;
    }
}// namespace astrokit
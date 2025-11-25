#include "ForceModel.h"

ForceModel::ForceModel(Planet& cb) : cb(cb), include_j2(true)
{
}

ForceModel::ForceModel(Planet& cb, bool include_j2) : cb(cb)
{
	set_include_j2(include_j2);
}

void ForceModel::set_include_j2(bool j2_included)
{
	this->include_j2 = j2_included;
}

Eigen::Vector<double, 6> ForceModel::eoms(double, const Eigen::Vector<double, 6>& state)
{
	Eigen::Vector<double, 6> dstate_dt_kep = astrokit::accel_kep(state, cb.get_mu());
	Eigen::Vector<double, 6> dstate_dt_j2 = astrokit::accel_j2(state, cb.get_mu(), cb.get_eq_radius(), cb.get_j2());
	//note: only want the accel components from accel_j2

	Eigen::Vector3d vel = dstate_dt_kep.segment<3>(0);
	Eigen::Vector3d accel = dstate_dt_kep.segment<3>(3);
	if (this->include_j2)
	{
		accel += dstate_dt_j2.segment<3>(3);
	}
	Eigen::Vector<double, 6> dstate_dt;
	dstate_dt << vel, accel;

	return dstate_dt;
}

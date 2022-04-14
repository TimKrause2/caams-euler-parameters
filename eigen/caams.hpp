#ifndef CAAMS_H
#define CAAMS_H

#include <Eigen/Dense>
#include "e2glm.h"

namespace Eigen
{
	typedef Matrix<double, 3, 4> Matrix3x4d;
}

namespace caams
{
	Eigen::Matrix3x4d G(Eigen::Vector4d const & p);
	Eigen::Matrix3x4d L(Eigen::Vector4d const & p);
	Eigen::Matrix3d SS(Eigen::Vector3d const & a);
	Eigen::Vector4d pAA(double angle, Eigen::Vector3d const & axis);
	Eigen::Matrix3d Ap(Eigen::Vector4d const & p);
	Eigen::Matrix4d a_plus(Eigen::Vector3d const & a);
	Eigen::Matrix4d a_minus(Eigen::Vector3d const & a);

    // some inertial tensor functions
	Eigen::Matrix3d J_p_cylinder_x_axis(double m, double r, double l);
	Eigen::Matrix3d J_p_cylinder_y_axis(double m, double r, double h);
	Eigen::Matrix3d J_p_cylinder_z_axis(double m, double r, double h);
	Eigen::Matrix3d J_p_sphere(double mass, double radius);
	Eigen::Matrix3d J_p_cuboid(double mass, double dx, double dy, double dz);

}

#endif

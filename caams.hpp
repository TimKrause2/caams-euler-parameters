#ifndef CAAMS_H
#define CAAMS_H

#include "matrix.hpp"

namespace caams
{
	matrix G(matrix const & p);
	matrix L(matrix const & p);
	matrix SS(matrix const & a);
	matrix pAA(double angle, matrix const & axis);
    matrix Ap(matrix const & p);
	matrix a_plus(matrix const & a);
	matrix a_minus(matrix const & a);
	double norm(matrix const & a);
	matrix x_axis(matrix const & a);
	matrix y_axis(matrix const & a);
	matrix z_axis(matrix const & a);
    matrix normalize(matrix const & a);
    long maxComponent(matrix const & a);

    // some inertial tensor functions
    matrix J_p_cylinder_x_axis(double m, double r, double l);
    matrix J_p_cylinder_y_axis(double m, double r, double h);
    matrix J_p_cylinder_z_axis(double m, double r, double h);
    matrix J_p_sphere(double mass, double radius);
    matrix J_p_cuboid(double mass, double dx, double dy, double dz);

}

#endif

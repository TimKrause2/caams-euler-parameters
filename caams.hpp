#include "matrix.hpp"

namespace caams
{
	matrix G(matrix const & p);
	matrix L(matrix const & p);
	matrix SS(matrix const & a);
	matrix pAA(double angle, matrix const & axis);
	matrix a_plus(matrix const & a);
	matrix a_minus(matrix const & a);
	double norm(matrix const & a);
	matrix x_axis(matrix const & a);
	matrix y_axis(matrix const & a);
	matrix z_axis(matrix const & a);
}
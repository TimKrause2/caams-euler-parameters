#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iostream>

int main(int argc, char **argv)
{
	Eigen::Vector3d v(1.0, 2.0, -3.0);
	int i_max;
	v.cwiseAbs().maxCoeff(&i_max);
	std::cout << "i_max:" << i_max << std::endl;
	return 0;

}

Eigen::MatrixXd some_func(void)
{
	Eigen::Matrix<double, 1, 7> r;
	return r;
}

#if 0

// This doesn't work. SparseMatrix doesn't support 'block'.
void test_func(void)
{
	Eigen::SparseMatrix<double> M(10,10);
	Eigen::Matrix<double,3,3> N;
	M.block<3,3>(0,0) = N;
}

#endif

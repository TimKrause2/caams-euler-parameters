#pragma once
#include <Eigen/Dense>
#include <glm/matrix.hpp>

template<typename T, int m, int n>
inline glm::mat<m, n, float, glm::precision::highp> E2GLM(const Eigen::Matrix<T, m, n>& em)
{
	glm::mat<m, n, float, glm::precision::highp> mat;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			mat[j][i] = em(i, j);
		}
	}
	return mat;
}

template<typename T, int m>
inline glm::vec<m, float, glm::precision::highp> E2GLM(const Eigen::Matrix<T, m, 1>& em)
{
	glm::vec<m, float, glm::precision::highp> v;
	for (int i = 0; i < m; ++i)
	{
		v[i] = em(i);
	}
	return v;
}

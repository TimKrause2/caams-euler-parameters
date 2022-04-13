#ifndef BODY_H
#define BODY_H

#include "caams.hpp"

class Body
{
public:
	Eigen::Vector3d r;
	Eigen::Vector4d p;
	Eigen::Vector3d r_dot;
	Eigen::Vector4d p_dot;
    double mass;
	Eigen::Matrix3d J_p;

	Eigen::Vector3d rk_r;
	Eigen::Vector4d rk_p;
	Eigen::Vector3d rk_r_dot;
	Eigen::Vector4d rk_p_dot;
	Eigen::Matrix3x4d k_r_dot;
	Eigen::Matrix4d k_p_dot;
	Eigen::Matrix3x4d k_r_ddot;
	Eigen::Matrix4d k_p_ddot;

    int eqn_index;
    int p_index;

    Body(
			Eigen::Vector3d r,
			Eigen::Vector4d p,
			Eigen::Vector3d r_dot,
			Eigen::Vector4d p_dot,
            double mass,
			Eigen::Matrix3d J_p);

    ~Body(void){}

	Eigen::Matrix3d N();
	Eigen::Matrix4d J_star();
	Eigen::Matrix<double,7,1> b_star();

    virtual void Draw(void)=0;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class DatumBody : public Body
{
public:
    DatumBody(void);
    ~DatumBody(void){}
    virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CylinderXaxis : public Body
{
    double radius;
    double length;
public:
    CylinderXaxis(
			Eigen::Vector3d r,
			Eigen::Vector4d p,
			Eigen::Vector3d r_dot,
			Eigen::Vector4d p_dot,
            double mass,
            double radius,
            double length);
    ~CylinderXaxis(void){}
    virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CylinderYaxis : public Body
{
    double radius;
    double length;
public:
    CylinderYaxis(
			Eigen::Vector3d r,
			Eigen::Vector4d p,
			Eigen::Vector3d r_dot,
			Eigen::Vector4d p_dot,
            double mass,
            double radius,
            double length);
    ~CylinderYaxis(void){}
    virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CylinderZaxis : public Body
{
    double radius;
    double length;
public:
    CylinderZaxis(
			Eigen::Vector3d r,
			Eigen::Vector4d p,
			Eigen::Vector3d r_dot,
			Eigen::Vector4d p_dot,
            double mass,
            double radius,
            double length);
    ~CylinderZaxis(void){}
    virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Cuboid : public Body
{
    double dx;
    double dy;
    double dz;
public:
    Cuboid(
			Eigen::Vector3d r,
			Eigen::Vector4d p,
			Eigen::Vector3d r_dot,
			Eigen::Vector4d p_dot,
            double mass,
            double dx,
            double dy,
            double dz);
    ~Cuboid(void){}
    virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


#endif

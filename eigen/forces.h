#ifndef FORCES_H
#define FORCES_H

#include "force_element.h"
#include "body.h"
#include "constraint.h"
#include "mbsystem.h"

class SystemGravityForce : public ForceElement
{
public:
	Eigen::Vector3d g;
    System &system;
	SystemGravityForce(Eigen::Vector3d g, System &system);
    ~SystemGravityForce(void){}
	void Apply(Eigen::VectorXd &y_rhs);
};

class BodyGravityForce : public ForceElement
{
public:
	Eigen::Vector3d g;
    Body *body;
	BodyGravityForce(Eigen::Vector3d g, Body *body);
    ~BodyGravityForce(void){}
	void Apply(Eigen::VectorXd &y_rhs);
};

class SystemNewtonianGravityForce : public ForceElement
{
public:
    double G_Newton;
    System &system;
    SystemNewtonianGravityForce(
            double G_Newton,
            System &system);
    ~SystemNewtonianGravityForce(){}
	void Apply(Eigen::VectorXd &y_rhs);
    void GravityForcePair(
            Body *body1,
            Body *body2,
			Eigen::VectorXd &y_rhs);
};



class LinearSpringDamperForce : public ForceElement
{
public:
    Body *body1;
    Body *body2;
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    double length;
    double k_spring;
    double k_damper;
    LinearSpringDamperForce(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p,
            double l0,
            double k_spring,
            double k_damper);
    ~LinearSpringDamperForce(){}
	void Apply(Eigen::VectorXd &y_rhs);
    void Draw(void);
};

class BodyLocalForce : public ForceElement
{
public:
    Body* body;
	Eigen::Vector3d f;
	BodyLocalForce(Body* body, Eigen::Vector3d f);
    ~BodyLocalForce(void){}
	void Apply(Eigen::VectorXd &y_rhs);
    void Draw(void);
};

class BodyGlobalForce : public ForceElement
{
public:
    Body* body;
	Eigen::Vector3d f;
	BodyGlobalForce(Body* body, Eigen::Vector3d f);
    ~BodyGlobalForce(void){}
	void Apply(Eigen::VectorXd &y_rhs);
    void Draw(void);
};

class BodyLocalTorque : public ForceElement
{
public:
    Body* body;
	Eigen::Vector3d n;
	BodyLocalTorque(Body *body, Eigen::Vector3d n);
    ~BodyLocalTorque(){}
	void Apply(Eigen::VectorXd &y_rhs);
    void Draw(void);
};

class BodyGlobalTorque : public ForceElement
{
public:
    Body* body;
	Eigen::Vector3d n;
	BodyGlobalTorque(Body *body, Eigen::Vector3d n);
    ~BodyGlobalTorque(){}
	void Apply(Eigen::VectorXd &y_rhs);
    void Draw(void);
};

class BodyDamping : public ForceElement
{
public:
    Body* body;
    Eigen::Matrix3d k_t;
    Eigen::Matrix3d k_r;
    BodyDamping(Body* body, Eigen::Matrix3d k_t, Eigen::Matrix3d k_r);
    ~BodyDamping(){}
	void Apply(Eigen::VectorXd &y_rhs);
    void Draw(void);
};


#endif

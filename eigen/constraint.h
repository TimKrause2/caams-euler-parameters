#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "caams.hpp"
#include "body.h"

class Constraint
{
public:
    Body *body1;
    Body *body2;
    int N_eqn;
    int eqn_index;
    Constraint(
            Body *body1,
            Body *body2);
	~Constraint(){}
	Eigen::VectorXd dPHI(void);
	Eigen::Vector3d h(const Eigen::Vector4d &p_dot,
					  const Eigen::Vector3d &s_p);
	virtual Eigen::MatrixXd Body1Jacobian(void)=0;
	virtual Eigen::MatrixXd Body2Jacobian(void)=0;
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void)=0;
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void)=0;
	virtual Eigen::VectorXd PHI(void)=0;
	virtual Eigen::VectorXd ModifiedGamma(void)=0;
    virtual void Draw(void)=0;
};

class SphericalJoint : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    SphericalJoint(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p);
	~SphericalJoint(){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};



class SphericalSphericalJoint : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    double length;
    SphericalSphericalJoint(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p,
            double length);
	~SphericalSphericalJoint(){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Normal1_1 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    Normal1_1(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p);
    ~Normal1_1(){};
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
   virtual void Draw(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Normal2_1 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s1B_p;
	Eigen::Vector3d s2B_p;
    Normal2_1(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s1B_p,
			Eigen::Vector3d s2B_p);
    ~Normal2_1(void){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
   virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Parallel1_2 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    Parallel1_2(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p);
    ~Parallel1_2(void){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Parallel2_2 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s1B_p;
	Eigen::Vector3d s2B_p;
    Parallel2_2(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s1B_p,
			Eigen::Vector3d s2B_p);
    ~Parallel2_2(void){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif

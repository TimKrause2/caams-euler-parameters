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
    ~Constraint(){};
    caams::matrix dPHI(void);
    caams::matrix h(const caams::matrix &p_dot, const caams::matrix &s_p);
    virtual caams::matrix Body1Jacobian(void)=0;
    virtual caams::matrix Body2Jacobian(void)=0;
    virtual caams::matrix Body1ModifiedJacobian(void)=0;
    virtual caams::matrix Body2ModifiedJacobian(void)=0;
    virtual caams::matrix PHI(void)=0;
    virtual caams::matrix ModifiedGamma(void)=0;
    virtual void Draw(void)=0;
};

class SphericalJoint : public Constraint
{
public:
    caams::matrix s1_p;
    caams::matrix s2_p;
    SphericalJoint(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s2_p);
    ~SphericalJoint(){};
    virtual caams::matrix Body1Jacobian(void);
    virtual caams::matrix Body2Jacobian(void);
    virtual caams::matrix Body1ModifiedJacobian(void);
    virtual caams::matrix Body2ModifiedJacobian(void);
    virtual caams::matrix PHI(void);
    virtual caams::matrix ModifiedGamma(void);
    virtual void Draw(void);

};



class SphericalSphericalJoint : public Constraint
{
public:
    caams::matrix s1_p;
    caams::matrix s2_p;
    double length;
    SphericalSphericalJoint(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s2_p,
            double length);
    ~SphericalSphericalJoint(){};
    virtual caams::matrix Body1Jacobian(void);
    virtual caams::matrix Body2Jacobian(void);
    virtual caams::matrix Body1ModifiedJacobian(void);
    virtual caams::matrix Body2ModifiedJacobian(void);
    virtual caams::matrix PHI(void);
    virtual caams::matrix ModifiedGamma(void);
    virtual void Draw(void);

};

class Normal1_1 : public Constraint
{
public:
    caams::matrix s1_p;
    caams::matrix s2_p;
    Normal1_1(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s2_p);
    ~Normal1_1(){};
    virtual caams::matrix Body1Jacobian(void);
    virtual caams::matrix Body2Jacobian(void);
    virtual caams::matrix Body1ModifiedJacobian(void);
    virtual caams::matrix Body2ModifiedJacobian(void);
    virtual caams::matrix PHI(void);
    virtual caams::matrix ModifiedGamma(void);
    virtual void Draw(void);

};

class Normal2_1 : public Constraint
{
public:
    caams::matrix s1_p;
    caams::matrix s1B_p;
    caams::matrix s2B_p;
    Normal2_1(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s1B_p,
            caams::matrix s2B_p);
    ~Normal2_1(void){}
    virtual caams::matrix Body1Jacobian(void);
    virtual caams::matrix Body2Jacobian(void);
    virtual caams::matrix Body1ModifiedJacobian(void);
    virtual caams::matrix Body2ModifiedJacobian(void);
    virtual caams::matrix PHI(void);
    virtual caams::matrix ModifiedGamma(void);
    virtual void Draw(void);
};

class Parallel1_2 : public Constraint
{
public:
    caams::matrix s1_p;
    caams::matrix s2_p;
    Parallel1_2(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s2_p);
    ~Parallel1_2(void){}
    virtual caams::matrix Body1Jacobian(void);
    virtual caams::matrix Body2Jacobian(void);
    virtual caams::matrix Body1ModifiedJacobian(void);
    virtual caams::matrix Body2ModifiedJacobian(void);
    virtual caams::matrix PHI(void);
    virtual caams::matrix ModifiedGamma(void);
    virtual void Draw(void);

};

class Parallel2_2 : public Constraint
{
public:
    caams::matrix s1_p;
    caams::matrix s1B_p;
    caams::matrix s2B_p;
    Parallel2_2(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s1B_p,
            caams::matrix s2B_p);
    ~Parallel2_2(void){}
    virtual caams::matrix Body1Jacobian(void);
    virtual caams::matrix Body2Jacobian(void);
    virtual caams::matrix Body1ModifiedJacobian(void);
    virtual caams::matrix Body2ModifiedJacobian(void);
    virtual caams::matrix PHI(void);
    virtual caams::matrix ModifiedGamma(void);
    virtual void Draw(void);
};






#endif




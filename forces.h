#ifndef FORCES_H
#define FORCES_H

#include "force_element.h"
#include "body.h"
#include "constraint.h"
#include "mbsystem.h"

class SystemGravityForce : public ForceElement
{
public:
    caams::matrix g;
    System &system;
    SystemGravityForce(caams::matrix g, System &system);
    ~SystemGravityForce(void){}
    void Apply(caams::matrix &y_rhs);
};

class BodyGravityForce : public ForceElement
{
public:
    caams::matrix g;
    Body *body;
    BodyGravityForce(caams::matrix g, Body *body);
    ~BodyGravityForce(void){}
    void Apply(caams::matrix &y_rhs);
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
    void Apply(caams::matrix &y_rhs);
    void GravityForcePair(
            Body *body1,
            Body *body2,
            caams::matrix &y_rhs);
};



class LinearSpringDamperForce : public ForceElement
{
public:
    Body *body1;
    Body *body2;
    caams::matrix s1_p;
    caams::matrix s2_p;
    double length;
    double k_spring;
    double k_damper;
    LinearSpringDamperForce(
            Body *body1,
            Body *body2,
            caams::matrix s1_p,
            caams::matrix s2_p,
            double l0,
            double k_spring,
            double k_damper);
    ~LinearSpringDamperForce(){}
    void Apply(caams::matrix &y_rhs);
    void Draw(void);
};



#endif

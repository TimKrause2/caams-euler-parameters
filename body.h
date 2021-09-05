#ifndef BODY_H
#define BODY_H

#include "caams.hpp"

class Body
{
public:
    caams::matrix r;
    caams::matrix p;
    caams::matrix r_dot;
    caams::matrix p_dot;
    double mass;
    caams::matrix J_p;

    caams::matrix rk_r;
    caams::matrix rk_p;
    caams::matrix rk_r_dot;
    caams::matrix rk_p_dot;
    caams::matrix k_r_dot;
    caams::matrix k_p_dot;
    caams::matrix k_r_ddot;
    caams::matrix k_p_ddot;

    int eqn_index;
    int p_index;

    Body(
            caams::matrix r,
            caams::matrix p,
            caams::matrix r_dot,
            caams::matrix p_dot,
            double mass,
            caams::matrix J_p);

    ~Body(void){}

    caams::matrix N();
    caams::matrix J_star();
    caams::matrix b_star();

    virtual void Draw(void)=0;
};

class DatumBody : public Body
{
public:
    DatumBody(void);
    ~DatumBody(void){}
    virtual void Draw(void);
};

class CylinderXaxis : public Body
{
    double radius;
    double length;
public:
    CylinderXaxis(
            caams::matrix r,
            caams::matrix p,
            caams::matrix r_dot,
            caams::matrix p_dot,
            double mass,
            double radius,
            double length);
    ~CylinderXaxis(void){}
    virtual void Draw(void);
};

class CylinderYaxis : public Body
{
    double radius;
    double length;
public:
    CylinderYaxis(
            caams::matrix r,
            caams::matrix p,
            caams::matrix r_dot,
            caams::matrix p_dot,
            double mass,
            double radius,
            double length);
    ~CylinderYaxis(void){}
    virtual void Draw(void);
};

class CylinderZaxis : public Body
{
    double radius;
    double length;
public:
    CylinderZaxis(
            caams::matrix r,
            caams::matrix p,
            caams::matrix r_dot,
            caams::matrix p_dot,
            double mass,
            double radius,
            double length);
    ~CylinderZaxis(void){}
    virtual void Draw(void);
};

class Cuboid : public Body
{
    double dx;
    double dy;
    double dz;
public:
    Cuboid(
            caams::matrix r,
            caams::matrix p,
            caams::matrix r_dot,
            caams::matrix p_dot,
            double mass,
            double dx,
            double dy,
            double dz);
    ~Cuboid(void){}
    virtual void Draw(void);
};










#endif

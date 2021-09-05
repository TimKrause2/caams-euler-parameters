#ifndef MBSYSTEM_H
#define MBSYSTEM_H

#include "caams.hpp"
#include "body.h"
#include "constraint.h"
#include "force_element.h"
#include <vector>

class SolverParameters
{
public:
    caams::matrix A;
    caams::matrix x;
    caams::matrix y_lhs;
    caams::matrix y_rhs;
    SolverParameters(long N_eqn):
        A(N_eqn,N_eqn,caams::zeros),
        x(N_eqn,1),
        y_lhs(N_eqn,1,caams::zeros),
        y_rhs(N_eqn,1,caams::zeros){};
    ~SolverParameters(void){}
};

class System
{
public:
    std::vector<Body*> bodies;
    std::vector<Constraint*> constraints;
    std::vector<ForceElement*> forces;

    long N_eqn;
    long p_index;
    long B_index;
    SolverParameters *sp;
    System(void);
    ~System(void);

    void AddBody(Body *body);
    void AddConstraint(Constraint *constraint);
    void AddForce(ForceElement *force);
    void InitializeSolver(void);

    void rkSolve(void);

    void rkPhase1State(double dt);
    void rkPhase2State(double dt);
    void rkPhase3State(double dt);
    void rkPhase4State(double dt);

    void rkPhase1Integrate(double dt);
    void rkPhase2Integrate(double dt);
    void rkPhase3Integrate(double dt);
    void rkPhase4Integrate(double dt);

    void Integrate(double dt);

    void Draw(void);
};

#endif

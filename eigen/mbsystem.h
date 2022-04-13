#ifndef MBSYSTEM_H
#define MBSYSTEM_H

#include "caams.hpp"
#include "body.h"
#include "constraint.h"
#include "force_element.h"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd EVector;
typedef Eigen::Triplet<double> T;

class SolverParameters
{
public:
	SpMat A;
	EVector x;
	EVector y_lhs;
	EVector y_rhs;
    SolverParameters(long N_eqn):
		A(N_eqn,N_eqn),
		x(N_eqn),
		y_lhs(N_eqn),
		y_rhs(N_eqn)
	{
		y_lhs = Eigen::VectorXd::Zero(N_eqn);
		y_rhs = Eigen::VectorXd::Zero(N_eqn);

	};
    ~SolverParameters(void){}
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
    long N_nzs;
    SolverParameters *sp;
    System(void);
    ~System(void);

    void AddBody(Body *body);
    void AddConstraint(Constraint *constraint);
    void AddForce(ForceElement *force);
    void InitializeSolver(void);
private:
    void rkSolve(void);

    void rkPhase1State(double dt);
    void rkPhase2State(double dt);
    void rkPhase3State(double dt);
    void rkPhase4State(double dt);

	void rkPhase1Integrate(void);
	void rkPhase2Integrate(void);
	void rkPhase3Integrate(void);
	void rkPhase4Integrate(void);

	void rkUpdateState(double dt);
public:
    void Integrate(double dt);

    void Draw(void);
};

#endif

#include "mbsystem.h"
#include <iostream>

#define STABILIZATION_ALPHA 8.0
#define STABILIZATION_BETA  8.0

System::System(void)
{
    sp = NULL;
}

System::~System(void)
{
    if(sp) delete sp;
}

void System::AddBody(Body *body)
{
    bodies.push_back(body);
}

void System::AddConstraint(Constraint *constraint)
{
    constraints.push_back(constraint);
}

void System::AddForce(ForceElement *force)
{
    forces.push_back(force);
}

void System::InitializeSolver(void)
{
	long eqn_index = 0;
    for(auto body:bodies){
        body->eqn_index = eqn_index;
        eqn_index += 7;
    }
    p_index = eqn_index;
    for(auto body:bodies){
        body->p_index = eqn_index;
        eqn_index++;
    }
    B_index = eqn_index;
    for(auto constraint:constraints){
        constraint->eqn_index = eqn_index;
        eqn_index += constraint->N_eqn;
    }

	N_eqn = eqn_index;

    std::cout << "Number of eqations:" << N_eqn << std::endl;

    sp = new SolverParameters(N_eqn);

    // estimate the number of non-zero entries
    int N_nzs = 0;
    N_nzs += bodies.size()*(3+16+2*4);

    for(auto constraint:constraints){
        if(constraint->body1->eqn_index!=0){
            N_nzs += 2*7*constraint->N_eqn;
        }
        if(constraint->body2->eqn_index!=0){
            N_nzs += 2*7*constraint->N_eqn;
        }
    }

    std::cout << "Number of non-zeros:" << N_nzs << std::endl;
}

void sub(const Eigen::MatrixXd &M, long row, long col, std::vector<T> &nzs){
	for(long i_col=0;i_col<M.cols();i_col++){
		for(long i_row=0;i_row<M.rows();i_row++){
			if(M(i_row,i_col) != 0.0){
				nzs.push_back(T(row+i_row,col+i_col,M(i_row,i_col)));
            }
        }
	}
}

void System::rkSolve(void)
{
    std::vector<T> nzs;
    nzs.reserve(N_nzs);

    // build the non-zero triplet list
    for(auto body:bodies){
        sub(body->N(),body->eqn_index,body->eqn_index,nzs);
        sub(body->J_star(),body->eqn_index+3,body->eqn_index+3,nzs);
		sub(body->rk_p.transpose(),body->p_index,body->eqn_index+3,nzs);
        sub(body->rk_p,body->eqn_index+3,body->p_index,nzs);
		sp->y_lhs.segment<7>(body->eqn_index) = body->b_star();
		sp->y_lhs.segment<1>(body->p_index) = body->rk_p_dot.transpose()*body->rk_p_dot;
		sp->y_rhs.segment<7>(body->eqn_index) = Eigen::Matrix<double,7,1>::Zero();
    }

    for(auto force:forces){
        force->Apply(sp->y_rhs);
    }

    for(auto constraint:constraints){
		if(constraint->body1->eqn_index>=0){
            sub(constraint->Body1ModifiedJacobian(),
                constraint->eqn_index,
                constraint->body1->eqn_index,
                nzs);
			sub(constraint->Body1Jacobian().transpose(),
                constraint->body1->eqn_index,
                constraint->eqn_index,
                nzs);
        }
		if(constraint->body2->eqn_index>=0){
            sub(constraint->Body2ModifiedJacobian(),
                constraint->eqn_index,
                constraint->body2->eqn_index,
                nzs);
			sub(constraint->Body2Jacobian().transpose(),
                constraint->body2->eqn_index,
                constraint->eqn_index,
                nzs);
        }
		Eigen::VectorXd gamma_t(constraint->ModifiedGamma());
		//std::cout << "gamma:" << std::endl << gamma_t << std::endl;
		gamma_t -= 2.0*STABILIZATION_ALPHA*constraint->dPHI();
		Eigen::VectorXd PHI(constraint->PHI());
        //PHI.print("PHI:");
		gamma_t -= STABILIZATION_BETA*STABILIZATION_BETA*PHI;
		sp->y_rhs.segment(constraint->eqn_index,constraint->N_eqn) = gamma_t;
    }

	sp->A.setFromTriplets(nzs.begin(), nzs.end());

    Eigen::SparseLU<SpMat> solver;

	solver.analyzePattern(sp->A);
	solver.factorize(sp->A);

	sp->x = solver.solve(sp->y_rhs - sp->y_lhs);

}

void System::rkPhase1State(double dt)
{
    for(auto body:bodies){
        body->rk_r = body->r;
        body->rk_p = body->p;
        body->rk_r_dot = body->r_dot;
        body->rk_p_dot = body->p_dot;
    }
}

void System::rkPhase2State(double dt)
{
    for(auto body:bodies){
		body->rk_r = body->r + dt/2.0*body->k_r_dot.col(0);
		body->rk_p = body->p + dt/2.0*body->k_p_dot.col(0);
		body->rk_p.normalize();
		body->rk_r_dot = body->r_dot + dt/2.0*body->k_r_ddot.col(0);
		body->rk_p_dot = body->p_dot + dt/2.0*body->k_p_ddot.col(0);
    }
}

void System::rkPhase3State(double dt)
{
    for(auto body:bodies){
		body->rk_r = body->r + dt/2.0*body->k_r_dot.col(1);
		body->rk_p = body->p + dt/2.0*body->k_p_dot.col(1);
		body->rk_p.normalize();
		body->rk_r_dot = body->r_dot + dt/2.0*body->k_r_ddot.col(1);
		body->rk_p_dot = body->p_dot + dt/2.0*body->k_p_ddot.col(1);
    }
}

void System::rkPhase4State(double dt)
{
    for(auto body:bodies){
		body->rk_r = body->r + dt*body->k_r_dot.col(2);
		body->rk_p = body->p + dt*body->k_p_dot.col(2);
		body->rk_p.normalize();
		body->rk_r_dot = body->r_dot + dt*body->k_r_ddot.col(2);
		body->rk_p_dot = body->p_dot + dt*body->k_p_ddot.col(2);
    }
}

void System::rkPhase1Integrate(void)
{
    for(auto body:bodies){
		body->k_r_dot.col(0) = body->rk_r_dot;
		body->k_p_dot.col(0) = body->rk_p_dot;
		body->k_r_ddot.col(0) = sp->x.segment<3>(body->eqn_index);
		body->k_p_ddot.col(0) = sp->x.segment<4>(body->eqn_index+3);
    }
}

void System::rkPhase2Integrate(void)
{
    for(auto body:bodies){
		body->k_r_dot.col(1) = body->rk_r_dot;
		body->k_p_dot.col(1) = body->rk_p_dot;
		body->k_r_ddot.col(1) = sp->x.segment<3>(body->eqn_index);
		body->k_p_ddot.col(1) = sp->x.segment<4>(body->eqn_index+3);
	}
}

void System::rkPhase3Integrate(void)
{
    for(auto body:bodies){
		body->k_r_dot.col(2) = body->rk_r_dot;
		body->k_p_dot.col(2) = body->rk_p_dot;
		body->k_r_ddot.col(2) = sp->x.segment<3>(body->eqn_index);
		body->k_p_ddot.col(2) = sp->x.segment<4>(body->eqn_index+3);
	}
}

void System::rkPhase4Integrate(void)
{
    for(auto body:bodies){
		body->k_r_dot.col(3) = body->rk_r_dot;
		body->k_p_dot.col(3) = body->rk_p_dot;
		body->k_r_ddot.col(3) = sp->x.segment<3>(body->eqn_index);
		body->k_p_ddot.col(3) = sp->x.segment<4>(body->eqn_index+3);
	}
}

void System::rkUpdateState(double dt)
{
	Eigen::Vector4d c(1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0);
	for(auto body:bodies){
		body->r += dt*(body->k_r_dot*c);
		body->p += dt*(body->k_p_dot*c);
		body->p.normalize();
		body->r_dot += dt*(body->k_r_ddot*c);
		body->p_dot += dt*(body->k_p_ddot*c);
		double sigma = (body->p_dot.transpose()*body->p)(0);
		body->p_dot -= sigma*body->p;
    }
}

void System::Integrate(double dt)
{
    if(!sp)return;

    rkPhase1State(dt);
    rkSolve();
	rkPhase1Integrate();

    rkPhase2State(dt);
    rkSolve();
	rkPhase2Integrate();

    rkPhase3State(dt);
    rkSolve();
	rkPhase3Integrate();

    rkPhase4State(dt);
    rkSolve();
	rkPhase4Integrate();

	rkUpdateState(dt);
}

void System::Draw(void)
{
    for(auto body:bodies){
        body->Draw();
    }
    for(auto constraint:constraints){
        constraint->Draw();
    }
    for(auto force:forces){
        force->Draw();
    }
}

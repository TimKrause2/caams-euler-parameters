#include "mbsystem.h"

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
    long eqn_index = 1;
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

    N_eqn = eqn_index-1;

    sp = new SolverParameters(N_eqn);

}

void System::rkSolve(void)
{
    for(auto body:bodies){
        sp->A.sub(body->N(),body->eqn_index,body->eqn_index);
        sp->A.sub(body->J_star(),body->eqn_index+3,body->eqn_index+3);
        sp->A.sub(~body->rk_p,body->p_index,body->eqn_index+3);
        sp->A.sub(body->rk_p,body->eqn_index+3,body->p_index);
        sp->y_lhs.sub(body->b_star(),body->eqn_index,1);
        sp->y_lhs.sub(~body->rk_p_dot*body->rk_p_dot,body->p_index,1);
        sp->y_rhs.sub(caams::matrix(3,1,0.0,0.0,0.0),body->eqn_index,1);
        sp->y_rhs.sub(caams::matrix(4,1,0.0,0.0,0.0,0.0),body->eqn_index+3,1);
    }

    for(auto force:forces){
        force->Apply(sp->y_rhs);
    }

    for(auto constraint:constraints){
        if(constraint->body1->eqn_index!=0){
            sp->A.sub(constraint->Body1ModifiedJacobian(),
                      constraint->eqn_index,
                      constraint->body1->eqn_index);
            sp->A.sub(~constraint->Body1Jacobian(),
                      constraint->body1->eqn_index,
                      constraint->eqn_index);
        }
        if(constraint->body2->eqn_index!=0){
            sp->A.sub(constraint->Body2ModifiedJacobian(),
                      constraint->eqn_index,
                      constraint->body2->eqn_index);
            sp->A.sub(~constraint->Body2Jacobian(),
                      constraint->body2->eqn_index,
                      constraint->eqn_index);
        }
        caams::matrix gamma_t(constraint->ModifiedGamma());
        //gamma_t.print("gamma:");
        gamma_t -= 2.0*STABILIZATION_ALPHA*constraint->dPHI();
        caams::matrix PHI(constraint->PHI());
        //PHI.print("PHI:");
        gamma_t -= STABILIZATION_BETA*STABILIZATION_BETA*PHI;
        sp->y_rhs.sub(gamma_t,constraint->eqn_index,1);
    }

    sp->x = sp->A.inverse()*(sp->y_rhs-sp->y_lhs);

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
        body->rk_r = body->r + dt/2.0*body->k_r_dot.sub(3,1,1,1);
        body->rk_p = normalize(body->p + dt/2.0*body->k_p_dot.sub(4,1,1,1));
        body->rk_r_dot = body->r_dot + dt/2.0*body->k_r_ddot.sub(3,1,1,1);
        body->rk_p_dot = body->p_dot + dt/2.0*body->k_p_ddot.sub(4,1,1,1);
    }
}

void System::rkPhase3State(double dt)
{
    for(auto body:bodies){
        body->rk_r = body->r + dt/2.0*body->k_r_dot.sub(3,1,1,2);
        body->rk_p = normalize(body->p + dt/2.0*body->k_p_dot.sub(4,1,1,2));
        body->rk_r_dot = body->r_dot + dt/2.0*body->k_r_ddot.sub(3,1,1,2);
        body->rk_p_dot = body->p_dot + dt/2.0*body->k_p_ddot.sub(4,1,1,2);
    }
}

void System::rkPhase4State(double dt)
{
    for(auto body:bodies){
        body->rk_r = body->r + dt*body->k_r_dot.sub(3,1,1,3);
        body->rk_p = normalize(body->p + dt*body->k_p_dot.sub(4,1,1,3));
        body->rk_r_dot = body->r_dot + dt*body->k_r_ddot.sub(3,1,1,3);
        body->rk_p_dot = body->p_dot + dt*body->k_p_ddot.sub(4,1,1,3);
    }
}

void System::rkPhase1Integrate(double dt)
{
    for(auto body:bodies){
        body->k_r_dot.sub(body->r_dot,1,1);
        body->k_p_dot.sub(body->p_dot,1,1);
        body->k_r_ddot.sub(sp->x.sub(3,1,body->eqn_index,1),1,1);
        body->k_p_ddot.sub(sp->x.sub(4,1,body->eqn_index+3,1),1,1);
    }
}

void System::rkPhase2Integrate(double dt)
{
    for(auto body:bodies){
        body->k_r_dot.sub(body->rk_r_dot,1,2);
        body->k_p_dot.sub(body->rk_p_dot,1,2);
        body->k_r_ddot.sub(sp->x.sub(3,1,body->eqn_index,1),1,2);
        body->k_p_ddot.sub(sp->x.sub(4,1,body->eqn_index+3,1),1,2);
    }
}

void System::rkPhase3Integrate(double dt)
{
    for(auto body:bodies){
        body->k_r_dot.sub(body->rk_r_dot,1,3);
        body->k_p_dot.sub(body->rk_p_dot,1,3);
        body->k_r_ddot.sub(sp->x.sub(3,1,body->eqn_index,1),1,3);
        body->k_p_ddot.sub(sp->x.sub(4,1,body->eqn_index+3,1),1,3);
    }
}

void System::rkPhase4Integrate(double dt)
{
    for(auto body:bodies){
        body->k_r_dot.sub(body->rk_r_dot,1,4);
        body->k_p_dot.sub(body->rk_p_dot,1,4);
        body->k_r_ddot.sub(sp->x.sub(3,1,body->eqn_index,1),1,4);
        body->k_p_ddot.sub(sp->x.sub(4,1,body->eqn_index+3,1),1,4);
        caams::matrix c(4,1,1.0,2.0,2.0,1.0);
        body->r += (dt/6.0)*(body->k_r_dot*c);
        body->p = caams::normalize(body->p + (dt/6.0)*(body->k_p_dot*c));
        body->r_dot += (dt/6.0)*(body->k_r_ddot*c);
        body->p_dot += (dt/6.0)*(body->k_p_ddot*c);
    }
}

void System::Integrate(double dt)
{
    if(!sp)return;

    rkPhase1State(dt);
    rkSolve();
    rkPhase1Integrate(dt);

    rkPhase2State(dt);
    rkSolve();
    rkPhase2Integrate(dt);

    rkPhase3State(dt);
    rkSolve();
    rkPhase3Integrate(dt);

    rkPhase4State(dt);
    rkSolve();
    rkPhase4Integrate(dt);
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

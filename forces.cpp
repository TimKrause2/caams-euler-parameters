#include "forces.h"
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <stdio.h>

SystemGravityForce::SystemGravityForce(
        caams::matrix g,
        System &system):
    g(g),
    system(system)
{

}

void SystemGravityForce::Apply(caams::matrix &y_rhs)
{
    for(auto body:system.bodies){
        caams::matrix f(body->mass*g);
        AccumulateForce(y_rhs,f,body->eqn_index);
    }
}

BodyGravityForce::BodyGravityForce(
        caams::matrix g,
        Body *body):
    g(g),
    body(body){

}

void BodyGravityForce::Apply(caams::matrix &y_rhs)
{
    caams::matrix f(body->mass*g);
    AccumulateForce(y_rhs, f, body->eqn_index);
}


SystemNewtonianGravityForce::SystemNewtonianGravityForce(
        double G_Newton,
        System &system):
    G_Newton(G_Newton),
    system(system)
{

}

void SystemNewtonianGravityForce::Apply(caams::matrix &y_rhs)
{
    long Nbodies = system.bodies.size();
    if(Nbodies<=1)return;
    for(long b1=0;b1<(Nbodies-1);b1++){
        for(long b2=b1+1;b2<Nbodies;b2++){
            GravityForcePair(
                        system.bodies[b1],
                        system.bodies[b2],
                        y_rhs);
        }
    }
}

void SystemNewtonianGravityForce::GravityForcePair(
        Body *body1,
        Body *body2,
        caams::matrix &y_rhs)
{
    caams::matrix r21(body2->rk_r - body1->rk_r);
    double mag_r21 = caams::norm(r21);
    caams::matrix r21n(1.0/mag_r21*r21);
    double mag_F1 = body1->mass*body2->mass*G_Newton/(mag_r21*mag_r21);
    caams::matrix F1(mag_F1*r21n);
    AccumulateForce(y_rhs,F1,body1->eqn_index);
    AccumulateForce(y_rhs,-F1,body2->eqn_index);

}



LinearSpringDamperForce::LinearSpringDamperForce(
        Body *body1,
        Body *body2,
        caams::matrix s1_p,
        caams::matrix s2_p,
        double length,
        double k_spring,
        double k_damper):
    body1(body1),
    body2(body2),
    s1_p(s1_p),
    s2_p(s2_p),
    length(length),
    k_spring(k_spring),
    k_damper(k_damper)
{

}

void LinearSpringDamperForce::Apply(caams::matrix &y_rhs)
{
    caams::matrix A1(caams::Ap(body1->rk_p));
    caams::matrix A2(caams::Ap(body2->rk_p));
    caams::matrix s1(A1*s1_p);
    caams::matrix s2(A2*s2_p);
    caams::matrix r1(body1->rk_r + s1);
    caams::matrix r2(body2->rk_r + s2);
    caams::matrix r21(r2-r1);
    double mag_r21 = caams::norm(r21);
    //printf("mag_r21:%lf\n",mag_r21);
    caams::matrix r21n(caams::normalize(r21));
    caams::matrix omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
    caams::matrix omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
    caams::matrix v1(body1->rk_r_dot + caams::SS(omega1)*s1);
    caams::matrix v2(body2->rk_r_dot + caams::SS(omega2)*s2);
    caams::matrix v21(v2-v1);
    caams::matrix v21_dot_r21n(~v21*r21n);
    double mag_F1 = (mag_r21-length)*k_spring + v21_dot_r21n.data[0]*k_damper;
    caams::matrix F1(mag_F1*r21n);
    caams::matrix F2(-F1);
    caams::matrix n1(caams::SS(s1)*F1);
    caams::matrix n2(caams::SS(s2)*F2);

    if(body1->eqn_index){
        AccumulateForce(y_rhs,F1,body1->eqn_index);
        AccumulateTorque(y_rhs,2.0*~caams::G(body1->rk_p)*n1,body1->eqn_index+3);
    }

    if(body2->eqn_index){
        AccumulateForce(y_rhs,F2,body2->eqn_index);
        AccumulateTorque(y_rhs,2.0*~caams::G(body2->rk_p)*n2,body2->eqn_index+3);
    }

}

void LinearSpringDamperForce::Draw(void)
{
    caams::matrix r1(body1->r + caams::Ap(body1->p)*s1_p);
    caams::matrix r2(body2->r + caams::Ap(body2->p)*s2_p);

    glBegin(GL_LINES);
        glColor3f(0.5f,0.5f,0.5f);
        glVertex3dv(body1->r.data);
        glVertex3dv(r1.data);
        glColor3f(0.75f,0.5f,0.0f);
        glVertex3dv(r1.data);
        glVertex3dv(r2.data);
        glColor3f(0.5f,0.5f,0.5f);
        glVertex3dv(r2.data);
        glVertex3dv(body2->r.data);
    glEnd();
}

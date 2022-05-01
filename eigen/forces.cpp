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
		Eigen::Vector3d g,
        System &system):
    g(g),
    system(system)
{

}

void SystemGravityForce::Apply(Eigen::VectorXd &y_rhs)
{
    for(auto body:system.bodies){
		Eigen::Vector3d f(body->mass*g);
        AccumulateForce(y_rhs,f,body->eqn_index);
    }
}

BodyGravityForce::BodyGravityForce(
		Eigen::Vector3d g,
        Body *body):
    g(g),
    body(body){

}

void BodyGravityForce::Apply(Eigen::VectorXd &y_rhs)
{
	Eigen::Vector3d f(body->mass*g);
    AccumulateForce(y_rhs, f, body->eqn_index);
}


SystemNewtonianGravityForce::SystemNewtonianGravityForce(
        double G_Newton,
        System &system):
    G_Newton(G_Newton),
    system(system)
{

}

void SystemNewtonianGravityForce::Apply(Eigen::VectorXd &y_rhs)
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
		Eigen::VectorXd &y_rhs)
{
	Eigen::Vector3d r21(body2->rk_r - body1->rk_r);
	double mag_r21 = r21.norm();
	Eigen::Vector3d r21n(1.0/mag_r21*r21);
    double mag_F1 = body1->mass*body2->mass*G_Newton/(mag_r21*mag_r21);
	Eigen::Vector3d F1(mag_F1*r21n);
    AccumulateForce(y_rhs,F1,body1->eqn_index);
    AccumulateForce(y_rhs,-F1,body2->eqn_index);

}



LinearSpringDamperForce::LinearSpringDamperForce(
        Body *body1,
        Body *body2,
		Eigen::Vector3d s1_p,
		Eigen::Vector3d s2_p,
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

void LinearSpringDamperForce::Apply(Eigen::VectorXd &y_rhs)
{
	Eigen::Matrix3d A1(caams::Ap(body1->rk_p));
	Eigen::Matrix3d A2(caams::Ap(body2->rk_p));
	Eigen::Vector3d s1(A1*s1_p);
	Eigen::Vector3d s2(A2*s2_p);
	Eigen::Vector3d r1(body1->rk_r + s1);
	Eigen::Vector3d r2(body2->rk_r + s2);
	Eigen::Vector3d r21(r2-r1);
	double mag_r21 = r21.norm();
    //printf("mag_r21:%lf\n",mag_r21);
	Eigen::Vector3d r21n = r21;
	r21n.normalize();
	Eigen::Vector3d omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
	Eigen::Vector3d omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
	Eigen::Vector3d v1(body1->rk_r_dot + caams::SS(omega1)*s1);
	Eigen::Vector3d v2(body2->rk_r_dot + caams::SS(omega2)*s2);
	Eigen::Vector3d v21(v2-v1);
	Eigen::Matrix<double,1,1> v21_dot_r21n = v21.transpose()*r21n;
	double mag_F1 = (mag_r21-length)*k_spring + v21_dot_r21n(0,0)*k_damper;
	Eigen::Vector3d F1(mag_F1*r21n);
	Eigen::Vector3d F2(-F1);
	Eigen::Vector3d n1(caams::SS(s1)*F1);
	Eigen::Vector3d n2(caams::SS(s2)*F2);

	if(body1->eqn_index>=0){
        AccumulateForce(y_rhs,F1,body1->eqn_index);
		AccumulateTorque(y_rhs,2.0*caams::G(body1->rk_p).transpose()*n1,body1->eqn_index+3);
    }

	if(body2->eqn_index>=0){
        AccumulateForce(y_rhs,F2,body2->eqn_index);
		AccumulateTorque(y_rhs,2.0*caams::G(body2->rk_p).transpose()*n2,body2->eqn_index+3);
    }

}

void LinearSpringDamperForce::Draw(void)
{
	Eigen::Vector3d r1(body1->r + caams::Ap(body1->p)*s1_p);
	Eigen::Vector3d r2(body2->r + caams::Ap(body2->p)*s2_p);

    glBegin(GL_LINES);
        glColor3f(0.5f,0.5f,0.5f);
		glVertex3dv(body1->r.data());
		glVertex3dv(r1.data());
        glColor3f(0.75f,0.5f,0.0f);
		glVertex3dv(r1.data());
		glVertex3dv(r2.data());
        glColor3f(0.5f,0.5f,0.5f);
		glVertex3dv(r2.data());
		glVertex3dv(body2->r.data());
    glEnd();
}

BodyLocalForce::BodyLocalForce(Body *body, Eigen::Vector3d f):
    body(body),
    f(f)
{

}

void BodyLocalForce::Apply(Eigen::VectorXd &y_rhs)
{
	Eigen::Vector3d f_g = caams::Ap(body->rk_p)*f;
    AccumulateForce(y_rhs, f_g, body->eqn_index);
}

void BodyLocalForce::Draw(void)
{

}

BodyGlobalForce::BodyGlobalForce(Body *body, Eigen::Vector3d f):
    body(body),
    f(f)
{

}

void BodyGlobalForce::Apply(Eigen::VectorXd &y_rhs)
{
    AccumulateForce(y_rhs, f, body->eqn_index);
}

void BodyGlobalForce::Draw(void)
{

}

BodyLocalTorque::BodyLocalTorque(Body *body, Eigen::Vector3d n):
    body(body),
    n(n)
{

}

void BodyLocalTorque::Apply(Eigen::VectorXd &y_rhs)
{
	AccumulateTorque(y_rhs,2.0*caams::L(body->rk_p).transpose()*n,body->eqn_index+3);
}

void BodyLocalTorque::Draw(void)
{

}

BodyGlobalTorque::BodyGlobalTorque(Body *body, Eigen::Vector3d n):
    body(body),
    n(n)
{

}

void BodyGlobalTorque::Apply(Eigen::VectorXd &y_rhs)
{
	AccumulateTorque(y_rhs,2.0*caams::G(body->rk_p).transpose()*n,body->eqn_index+3);
}

void BodyGlobalTorque::Draw(void)
{

}

BodyDamping::BodyDamping(Body *body, Eigen::Matrix3d k_t, Eigen::Matrix3d k_r):
    body(body),
    k_t(k_t),
    k_r(k_r)
{

}

void BodyDamping::Apply(Eigen::VectorXd &y_rhs)
{
    Eigen::Vector3d omega_p(2.0*caams::L(body->rk_p)*body->rk_p_dot);
    Eigen::Vector3d n_p = -k_r * omega_p;
    AccumulateTorque(y_rhs,2.0*caams::L(body->rk_p).transpose()*n_p,body->eqn_index+3);
    Eigen::Matrix3d A(caams::Ap(body->rk_p));
    Eigen::Vector3d r_dot_p = A.transpose()*body->rk_r_dot;
    Eigen::Vector3d F_p = -k_t * r_dot_p;
    Eigen::Vector3d F = A*F_p;
    AccumulateForce(y_rhs, F, body->eqn_index);
}

void BodyDamping::Draw(void)
{

}

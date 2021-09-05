#include "constraint.h"
#include <stdio.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>

Constraint::Constraint(
        Body *body1,
        Body *body2):
    body1(body1),
    body2(body2)
{

}

caams::matrix Constraint::dPHI(void)
{
    caams::matrix q1_dot(7,1);
    caams::matrix q2_dot(7,1);

    q1_dot.sub(body1->rk_r_dot,1,1);
    q1_dot.sub(body1->rk_p_dot,4,1);
    q2_dot.sub(body2->rk_r_dot,1,1);
    q2_dot.sub(body2->rk_p_dot,4,1);

    return Body1Jacobian()*q1_dot + Body2Jacobian()*q2_dot;
}

caams::matrix Constraint::h(const caams::matrix &p_dot, const caams::matrix &s_p)
{
    return -2.0*caams::G(p_dot)*~caams::L(p_dot)*s_p;
}

SphericalJoint::SphericalJoint(
        Body *body1,
        Body *body2,
        caams::matrix s1_p,
        caams::matrix s2_p):
    Constraint(body1,body2),
    s1_p(s1_p),
    s2_p(s2_p)
{
    N_eqn = 3;
}

caams::matrix SphericalJoint::Body1Jacobian(void)
{
    caams::matrix r(3,7);
    r.sub(caams::matrix(3,3,caams::init_identity),1,1);
    r.sub(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*~body1->rk_p),1,4);
    return r;
}

caams::matrix SphericalJoint::Body2Jacobian(void)
{
    caams::matrix r(3,7);
    r.sub(-caams::matrix(3,3,caams::init_identity),1,1);
    r.sub(-2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)+s2_p*~body2->rk_p),1,4);
    return r;
}

caams::matrix SphericalJoint::Body1ModifiedJacobian(void)
{
    caams::matrix r(3,7);
    r.sub(caams::matrix(3,3,caams::init_identity),1,1);
    r.sub(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)),1,4);
    return r;
}

caams::matrix SphericalJoint::Body2ModifiedJacobian(void)
{
    caams::matrix r(3,7);
    r.sub(-caams::matrix(3,3,caams::init_identity),1,1);
    r.sub(-2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)),1,4);
    return r;
}

caams::matrix SphericalJoint::PHI(void)
{
    return body1->rk_r + caams::Ap(body1->rk_p)*s1_p
           - body2->rk_r - caams::Ap(body2->rk_p)*s2_p;
}

caams::matrix SphericalJoint::ModifiedGamma(void)
{
    return h(body1->rk_p_dot,s1_p) - h(body2->rk_p_dot,s2_p);
}

void SphericalJoint::Draw(void)
{
    caams::matrix s1(body1->r + caams::Ap(body1->p)*s1_p);
    caams::matrix s2(body2->r + caams::Ap(body2->p)*s2_p);

    glBegin(GL_LINES);
        glColor3f(1.0f,1.0f,0.0f);
        glVertex3dv(body1->r.data);
        glVertex3dv(s1.data);
        glColor3f(0.0f,1.0f,1.0f);
        glVertex3dv(body2->r.data);
        glVertex3dv(s2.data);
    glEnd();
}

SphericalSphericalJoint::SphericalSphericalJoint(
        Body *body1,
        Body *body2,
        caams::matrix s1_p,
        caams::matrix s2_p,
        double length):
    Constraint(body1,body2),
    s1_p(s1_p),
    s2_p(s2_p),
    length(length)
{
    N_eqn = 1;
}

caams::matrix SphericalSphericalJoint::Body1Jacobian(void)
{
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(1,7);
    r.sub(-2.0*~d,1,1);
    caams::matrix B1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)
                     + s1_p*~body1->rk_p));
    r.sub(-2.0*~d*B1,1,4);
    return r;
}

caams::matrix SphericalSphericalJoint::Body2Jacobian(void)
{
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(1,7);
    r.sub(2.0*~d,1,1);
    caams::matrix B2(2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)
                     + s2_p*~body2->rk_p));
    r.sub(2.0*~d*B2,1,4);
    return r;
}

caams::matrix SphericalSphericalJoint::Body1ModifiedJacobian(void)
{
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(1,7);
    r.sub(-2.0*~d,1,1);
    r.sub(-4.0*~d*caams::G(body1->rk_p)*caams::a_minus(s1_p),1,4);
    return r;
}

caams::matrix SphericalSphericalJoint::Body2ModifiedJacobian(void)
{
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(1,7);
    r.sub(2.0*~d,1,1);
    r.sub(4.0*~d*caams::G(body2->rk_p)*caams::a_minus(s2_p),1,4);
    return r;
}

caams::matrix SphericalSphericalJoint::PHI(void)
{
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1_p);

    return ~d*d - caams::matrix(1,1,length*length);
}

caams::matrix SphericalSphericalJoint::ModifiedGamma(void)
{
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1_p);

    caams::matrix omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
    caams::matrix omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);

    caams::matrix d_dot(
                body2->rk_r_dot + caams::SS(omega2)*caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r_dot - caams::SS(omega1)*caams::Ap(body1->rk_p)*s1_p);

    return 2.0*(~d*(h(body2->rk_p_dot,s2_p)-h(body1->rk_p_dot,s1_p))-~d_dot*d_dot);
}

void SphericalSphericalJoint::Draw(void)
{
    caams::matrix s1(body1->r + caams::Ap(body1->p)*s1_p);
    caams::matrix s2(body2->r + caams::Ap(body2->p)*s2_p);

    caams::matrix d(s2-s1);
    double l = caams::norm(d);
    printf("SpericalSphericalJoint length:%lf\n",l);

    glBegin(GL_LINES);
        glColor3f(1.0f,0.5f,0.5f);
        glVertex3dv(body1->r.data);
        glVertex3dv(s1.data);
        glColor3f(0.5f,0.5f,0.5f);
        glVertex3dv(s1.data);
        glVertex3dv(s2.data);
        glColor3f(0.5f,1.0f,0.5f);
        glVertex3dv(s2.data);
        glVertex3dv(body2->r.data);
    glEnd();
}

Normal1_1::Normal1_1(
                Body *body1,
                Body *body2,
                caams::matrix s1_p,
                caams::matrix s2_p):
    Constraint(body1,body2),
    s1_p(s1_p),
    s2_p(s2_p)
{
    N_eqn = 1;
}

caams::matrix Normal1_1::Body1Jacobian(void)
{
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix C1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*~body1->rk_p));
    caams::matrix r(1,7,caams::zeros);
    r.sub(~s2*C1,1,4);
    return r;
}

caams::matrix Normal1_1::Body2Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix C2(2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)+s2_p*~body2->rk_p));
    caams::matrix r(1,7,caams::zeros);
    r.sub(~s1*C2,1,4);
    return r;
}

caams::matrix Normal1_1::Body1ModifiedJacobian(void)
{
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix r(1,7,caams::zeros);
    r.sub(~s2*caams::G(body1->rk_p)*caams::a_minus(s1_p),1,4);
    return r;
}

caams::matrix Normal1_1::Body2ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(1,7,caams::zeros);
    r.sub(~s1*caams::G(body2->rk_p)*caams::a_minus(s2_p),1,4);
    return r;
}

caams::matrix Normal1_1::PHI(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    return ~s1*s2;
}

caams::matrix Normal1_1::ModifiedGamma(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
    caams::matrix omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
    caams::matrix s1_dot(caams::SS(omega1)*s1);
    caams::matrix s2_dot(caams::SS(omega2)*s2);
    return ~s1*h(body2->rk_p_dot,s2_p) + ~s2*h(body1->rk_p_dot,s1_p)
            - 2.0*~s1_dot*s2_dot;
}


void Normal1_1::Draw(void)
{

}

Normal2_1::Normal2_1(
        Body *body1,
        Body *body2,
        caams::matrix s1_p,
        caams::matrix s1B_p,
        caams::matrix s2B_p):
    Constraint(body1,body2),
    s1_p(s1_p),
    s1B_p(s1B_p),
    s2B_p(s2B_p)
{
    N_eqn = 1;
}

caams::matrix Normal2_1::Body1Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(body2->rk_r + caams::Ap(body2->rk_p)*s2B_p
                    -body1->rk_r - caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix B1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1B_p)+s1B_p*~body1->rk_p));
    caams::matrix C1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*~body1->rk_p));
    caams::matrix r(1,7);
    r.sub(-~s1,1,1);
    r.sub(-~s1*B1 + ~d*C1,1,4);
    return r;
}

caams::matrix Normal2_1::Body2Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix B2(2.0*(caams::G(body2->rk_p)*caams::a_minus(s2B_p)+s2B_p*~body1->rk_p));
    caams::matrix r(1,7);
    r.sub(~s1,1,1);
    r.sub(~s1*B2,1,4);
    return r;
}

caams::matrix Normal2_1::Body1ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(body2->rk_r + caams::Ap(body2->rk_p)*s2B_p
                    -body1->rk_r - caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix r(1,7);
    r.sub(-~s1,1,1);
    r.sub(-2.0*~s1*caams::G(body1->rk_p)*caams::a_minus(s1B_p)
          + 2.0*~d*caams::G(body1->rk_p)*caams::a_minus(s1_p),1,4);
    return r;
}

caams::matrix Normal2_1::Body2ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(1,7);
    r.sub(~s1,1,1);
    r.sub(2.0*~s1*caams::G(body2->rk_p)*caams::a_minus(s2B_p),1,4);
    return r;
}

caams::matrix Normal2_1::PHI(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(body2->rk_r + caams::Ap(body2->rk_p)*s2B_p
                    -body1->rk_r - caams::Ap(body1->rk_p)*s1B_p);
    return ~s1*d;
}

caams::matrix Normal2_1::ModifiedGamma(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(body2->rk_r + caams::Ap(body2->rk_p)*s2B_p
                    -body1->rk_r - caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
    caams::matrix omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
    caams::matrix s1_dot(caams::SS(omega1)*s1);
    caams::matrix d_dot(
                body2->rk_r_dot + caams::SS(omega2)*caams::Ap(body2->rk_p)*s2B_p
                - body1->rk_r_dot - caams::SS(omega1)*caams::Ap(body1->rk_p)*s1B_p);

    return ~s1*(h(body2->rk_p_dot,s2B_p)-h(body1->rk_p_dot,s1B_p))
            + ~d*h(body1->rk_p_dot,s1_p)
            - 2.0*~s1_dot*d_dot;
}

void Normal2_1::Draw(void)
{

}

Parallel1_2::Parallel1_2(
        Body *body1,
        Body *body2,
        caams::matrix s1_p,
        caams::matrix s2_p):
    Constraint(body1,body2),
    s1_p(s1_p),
    s2_p(s2_p)
{
    N_eqn = 2;
}

caams::matrix Parallel1_2::Body1Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix C1(caams::G(body1->rk_p)*caams::a_minus(s1_p) + s1_p*~body1->rk_p);
    caams::matrix r(2,7,caams::zeros);
    caams::matrix PHI_p(-caams::SS(s2)*C1);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_p.sub(2,4,2,1),1,4);
        break;
    case 2:
        r.sub(PHI_p.sub(1,4,1,1),1,4);
        r.sub(PHI_p.sub(1,4,3,1),2,4);
        break;
    case 3:
        r.sub(PHI_p.sub(2,4,1,1),1,4);
        break;
    }
    return r;
}

caams::matrix Parallel1_2::Body2Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix C2(caams::G(body2->rk_p)*caams::a_minus(s2_p) + s2_p*~body2->rk_p);
    caams::matrix r(2,7,caams::zeros);
    caams::matrix PHI_p(caams::SS(s1)*C2);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_p.sub(2,4,2,1),1,4);
        break;
    case 2:
        r.sub(PHI_p.sub(1,4,1,1),1,4);
        r.sub(PHI_p.sub(1,4,3,1),2,4);
        break;
    case 3:
        r.sub(PHI_p.sub(2,4,1,1),1,4);
        break;
    }
    return r;
}

caams::matrix Parallel1_2::Body1ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix r(2,7,caams::zeros);
    caams::matrix PHI_p(-2.0*caams::SS(s2)*caams::G(body1->rk_p)*caams::a_minus(s1_p));
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_p.sub(2,4,2,1),1,4);
        break;
    case 2:
        r.sub(PHI_p.sub(1,4,1,1),1,4);
        r.sub(PHI_p.sub(1,4,3,1),2,4);
        break;
    case 3:
        r.sub(PHI_p.sub(2,4,1,1),1,4);
        break;
    }
    return r;
}

caams::matrix Parallel1_2::Body2ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(2,7,caams::zeros);
    caams::matrix PHI_p(2.0*caams::SS(s1)*caams::G(body2->rk_p)*caams::a_minus(s2_p));
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_p.sub(2,4,2,1),1,4);
        break;
    case 2:
        r.sub(PHI_p.sub(1,4,1,1),1,4);
        r.sub(PHI_p.sub(1,4,3,1),2,4);
        break;
    case 3:
        r.sub(PHI_p.sub(2,4,1,1),1,4);
        break;
    }
    return r;
}

caams::matrix Parallel1_2::PHI(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix r_all(caams::SS(s1)*s2);
    caams::matrix r(2,1);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(r_all.sub(2,1,2,1),1,1);
        break;
    case 2:
        r.sub(r_all.sub(1,1,1,1),1,1);
        r.sub(r_all.sub(1,1,3,1),2,1);
        break;
    case 3:
        r.sub(r_all.sub(2,1,1,1),1,1);
        break;
    }
    return r;
}

caams::matrix Parallel1_2::ModifiedGamma(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix s2(caams::Ap(body2->rk_p)*s2_p);
    caams::matrix omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
    caams::matrix omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
    caams::matrix s1_dot(caams::SS(omega1)*s1);
    caams::matrix s2_dot(caams::SS(omega2)*s2);
    caams::matrix r_all(caams::SS(s1)*h(body2->rk_p_dot,s2_p)
                        -caams::SS(s2)*h(body1->rk_p_dot,s1_p)
                        -2.0*caams::SS(s1_dot)*s2_dot);
    caams::matrix r(2,1);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(r_all.sub(2,1,2,1),1,1);
        break;
    case 2:
        r.sub(r_all.sub(1,1,1,1),1,1);
        r.sub(r_all.sub(1,1,3,1),2,1);
        break;
    case 3:
        r.sub(r_all.sub(2,1,1,1),1,1);
        break;
    }
    return r;
}

void Parallel1_2::Draw(void)
{

}

Parallel2_2::Parallel2_2(
        Body *body1,
        Body *body2,
        caams::matrix s1_p,
        caams::matrix s1B_p,
        caams::matrix s2B_p):
    Constraint(body1,body2),
    s1_p(s1_p),
    s1B_p(s1B_p),
    s2B_p(s2B_p)
{
    N_eqn = 2;
}

caams::matrix Parallel2_2::Body1Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix r(2,7);
    caams::matrix PHI_q(3,7);
    PHI_q.sub(-caams::SS(s1),1,1);
    caams::matrix B1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1B_p) + s1B_p*~body1->rk_p));
    caams::matrix C1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p) + s1_p*~body1->rk_p));
    PHI_q.sub(-caams::SS(s1)*B1 - caams::SS(d)*C1,1,4);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_q.sub(2,7,2,1),1,1);
        break;
    case 2:
        r.sub(PHI_q.sub(1,7,1,1),1,1);
        r.sub(PHI_q.sub(1,7,3,1),2,1);
        break;
    case 3:
        r.sub(PHI_q.sub(2,7,1,1),1,1);
        break;
    }
    return r;
}

caams::matrix Parallel2_2::Body2Jacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(2,7);
    caams::matrix PHI_q(3,7);
    PHI_q.sub(caams::SS(s1),1,1);
    caams::matrix B2(2.0*(caams::G(body2->rk_p)*caams::a_minus(s2B_p) + s2B_p*~body2->rk_p));
    PHI_q.sub(caams::SS(s1)*B2,1,4);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_q.sub(2,7,2,1),1,1);
        break;
    case 2:
        r.sub(PHI_q.sub(1,7,1,1),1,1);
        r.sub(PHI_q.sub(1,7,3,1),2,1);
        break;
    case 3:
        r.sub(PHI_q.sub(2,7,1,1),1,1);
        break;
    }
    return r;
}

caams::matrix Parallel2_2::Body1ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix r(2,7);
    caams::matrix PHI_q(3,7);
    PHI_q.sub(-caams::SS(s1),1,1);
    PHI_q.sub(-2.0*caams::SS(s1)*caams::G(body1->rk_p)*caams::a_minus(s1B_p)
              -2.0*caams::SS(d)*caams::G(body1->rk_p)*caams::a_minus(s1_p),1,4);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_q.sub(2,7,2,1),1,1);
        break;
    case 2:
        r.sub(PHI_q.sub(1,7,1,1),1,1);
        r.sub(PHI_q.sub(1,7,3,1),2,1);
        break;
    case 3:
        r.sub(PHI_q.sub(2,7,1,1),1,1);
        break;
    }
    return r;
}

caams::matrix Parallel2_2::Body2ModifiedJacobian(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix r(2,7);
    caams::matrix PHI_q(3,7);
    PHI_q.sub(caams::SS(s1),1,1);
    PHI_q.sub(2.0*caams::SS(s1)*caams::G(body2->rk_p)*caams::a_minus(s2B_p),1,4);
    switch(caams::maxComponent(s1)){
    case 1:
        r.sub(PHI_q.sub(2,7,2,1),1,1);
        break;
    case 2:
        r.sub(PHI_q.sub(1,7,1,1),1,1);
        r.sub(PHI_q.sub(1,7,3,1),2,1);
        break;
    case 3:
        r.sub(PHI_q.sub(2,7,1,1),1,1);
        break;
    }
    return r;
}

caams::matrix Parallel2_2::PHI(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix r_all(caams::SS(s1)*d);
    caams::matrix r(2,1);
    switch(caams::maxComponent(s1)){
    case 1:
        r.data[0] = r_all.data[1];
        r.data[1] = r_all.data[2];
        break;
    case 2:
        r.data[0] = r_all.data[0];
        r.data[1] = r_all.data[2];
        break;
    case 3:
        r.data[0] = r_all.data[0];
        r.data[1] = r_all.data[1];
        break;
    }
    return r;
}

caams::matrix Parallel2_2::ModifiedGamma(void)
{
    caams::matrix s1(caams::Ap(body1->rk_p)*s1_p);
    caams::matrix omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
    caams::matrix s1_dot(caams::SS(omega1)*s1);
    caams::matrix d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
    caams::matrix omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
    caams::matrix d_dot(
                body2->rk_r_dot + caams::SS(omega2)*caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r_dot - caams::SS(omega1)*caams::Ap(body1->rk_p)*s1B_p);

    caams::matrix r_all(caams::SS(s1)*(h(body2->rk_p_dot,s2B_p)-h(body1->rk_p_dot,s1B_p))
                        -caams::SS(d)*h(body1->rk_p_dot,s1_p)
                        -2.0*caams::SS(s1_dot)*d_dot);
    caams::matrix r(2,1);
    switch(caams::maxComponent(s1)){
    case 1:
        r.data[0] = r_all.data[1];
        r.data[1] = r_all.data[2];
        break;
    case 2:
        r.data[0] = r_all.data[0];
        r.data[1] = r_all.data[2];
        break;
    case 3:
        r.data[0] = r_all.data[0];
        r.data[1] = r_all.data[1];
        break;
    }
    return r;
}

void Parallel2_2::Draw(void)
{

}

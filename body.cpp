#include "body.h"
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>

Body::Body(
        caams::matrix r,
        caams::matrix p,
        caams::matrix r_dot,
        caams::matrix p_dot,
        double mass,
        caams::matrix J_p):
    r(r),
    p(p),
    r_dot(r_dot),
    p_dot(p_dot),
    mass(mass),
    J_p(J_p),
    rk_r(3,1,0.0,0.0,0.0),
    rk_p(4,1,1.0,0.0,0.0,0.0),
    rk_r_dot(3,1,0.0,0.0,0.0),
    rk_p_dot(4,1,0.0,0.0,0.0,0.0),
    k_r_dot(3,4),
    k_p_dot(4,4),
    k_r_ddot(3,4),
    k_p_ddot(4,4)
{
    eqn_index = 0;
}

caams::matrix Body::N(void)
{
    return mass*caams::matrix(3,3,caams::init_identity);
}

caams::matrix Body::J_star(void)
{
    return 4.0*~caams::L(rk_p)*J_p*caams::L(rk_p);
}

caams::matrix Body::b_star(void)
{
    caams::matrix r(7,1,caams::zeros);
    r.sub(8.0*~caams::L(rk_p_dot)*J_p*caams::L(rk_p)*rk_p_dot,4,1);
    return r;
}

DatumBody::DatumBody(void):
    Body(
        caams::matrix(3,1,0.0,0.0,0.0),
        caams::matrix(4,1,1.0,0.0,0.0,0.0),
        caams::matrix(3,1,0.0,0.0,0.0),
        caams::matrix(4,1,0.0,0.0,0.0,0.0),
        1.0,
        caams::matrix(3,3,caams::init_identity))
{

}

void DatumBody::Draw(void)
{
    caams::matrix A(caams::Ap(p));
    glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
    glm::dmat4x4 Abody(Aglm);
    glm::dvec3 rglm = glm::make_vec3(r.data);
    glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dmat4x4 Amodel = Atrans*Abody;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
    glBegin(GL_LINES);
        glColor3f(1.0f,0.0f,0.0f);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(1.0f,0.0f,0.0f);
        glColor3f(0.0f,1.0f,0.0f);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(0.0f,1.0f,0.0f);
        glColor3f(0.0f,0.0f,1.0f);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(0.0f,0.0f,1.0f);
    glEnd();
    glPopMatrix();

}




CylinderXaxis::CylinderXaxis(
        caams::matrix r,
        caams::matrix p,
        caams::matrix r_dot,
        caams::matrix p_dot,
        double mass,
        double radius,
        double length):
    Body(r,p,r_dot,p_dot,mass,caams::J_p_cylinder_x_axis(mass,radius,length)),
    radius(radius),
    length(length)
{

}

void CylinderXaxis::Draw(void)
{
    caams::matrix A(caams::Ap(p));
    glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
    glm::dmat4x4 Abody(Aglm);
    glm::dvec3 rglm = glm::make_vec3(r.data);
    glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 radj(0.0,0.0,-length/2.0);
    glm::dmat4x4 Aadj = glm::rotate(M_PI/2.0,glm::dvec3(0.0,1.0,0.0))*
                        glm::translate(radj);
    glm::dmat4x4 Amodel = Atrans*Abody*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
    glutWireCylinder(radius,length,(GLint)5,(GLint)1);
    glPopMatrix();
}


CylinderYaxis::CylinderYaxis(
        caams::matrix r,
        caams::matrix p,
        caams::matrix r_dot,
        caams::matrix p_dot,
        double mass,
        double radius,
        double length):
    Body(r,p,r_dot,p_dot,mass,caams::J_p_cylinder_y_axis(mass,radius,length)),
    radius(radius),
    length(length)
{

}

void CylinderYaxis::Draw(void)
{
    caams::matrix A(caams::Ap(p));
    glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
    glm::dmat4x4 Abody(Aglm);
    glm::dvec3 rglm = glm::make_vec3(r.data);
    glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 radj(0.0,0.0,-length/2.0);
    glm::dmat4x4 Aadj = glm::rotate(-M_PI/2.0,glm::dvec3(1.0,0.0,0.0))*
                        glm::translate(radj);
    glm::dmat4x4 Amodel = Atrans*Abody*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
    glutWireCylinder(radius,length,(GLint)5,(GLint)1);
    glPopMatrix();
}

CylinderZaxis::CylinderZaxis(
        caams::matrix r,
        caams::matrix p,
        caams::matrix r_dot,
        caams::matrix p_dot,
        double mass,
        double radius,
        double length):
    Body(r,p,r_dot,p_dot,mass,caams::J_p_cylinder_z_axis(mass,radius,length)),
    radius(radius),
    length(length)
{

}

void CylinderZaxis::Draw(void)
{
    caams::matrix A(caams::Ap(p));
    glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
    glm::dmat4x4 Abody(Aglm);
    glm::dvec3 rglm = glm::make_vec3(r.data);
    glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 radj(0.0,0.0,-length/2.0);
    glm::dmat4x4 Aadj = glm::translate(radj);
    glm::dmat4x4 Amodel = Atrans*Abody*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
    glutWireCylinder(radius,length,(GLint)5,(GLint)1);
    glPopMatrix();
}

Cuboid::Cuboid(
        caams::matrix r,
        caams::matrix p,
        caams::matrix r_dot,
        caams::matrix p_dot,
        double mass,
        double dx,
        double dy,
        double dz):
    Body(r,p,r_dot,p_dot,mass,caams::J_p_cuboid(mass,dx,dy,dz)),
    dx(dx),
    dy(dy),
    dz(dz)
{

}

void Cuboid::Draw(void)
{
    caams::matrix A(caams::Ap(p));
    glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
    glm::dmat4x4 Abody(Aglm);
    glm::dvec3 rglm = glm::make_vec3(r.data);
    glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 s(dx,dy,dz);
    glm::dmat4x4 Ascale = glm::scale(s);
    glm::dmat4x4 Amodel = Atrans*Abody*Ascale;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
    glutWireCube(1.0);
    glPopMatrix();

}



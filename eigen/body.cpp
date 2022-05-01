#include "body.h"
#include "primitives.h"
#include <iostream>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>

Body::Body(
		Eigen::Vector3d r,
		Eigen::Vector4d p,
		Eigen::Vector3d r_dot,
		Eigen::Vector4d p_dot,
		double mass,
		Eigen::Matrix3d J_p):
    r(r),
    p(p),
    r_dot(r_dot),
    p_dot(p_dot),
    mass(mass),
	J_p(J_p),
	rk_r(r),
	rk_p(p),
	rk_r_dot(r_dot),
	rk_p_dot(p_dot)
{
	eqn_index = -1;
}

Eigen::Matrix3d Body::N(void)
{
	return mass*Eigen::Matrix3d::Identity();
}

Eigen::Matrix4d Body::J_star(void)
{
	return 4.0*caams::L(rk_p).transpose()*J_p*caams::L(rk_p);
}

Eigen::Matrix<double,7,1> Body::b_star(void)
{
	Eigen::Matrix<double,7,1> r;
	r.head<3>() = Eigen::Vector3d::Zero();
	r.tail<4>() = 8.0*caams::L(rk_p_dot).transpose()*J_p*caams::L(rk_p)*rk_p_dot;
    return r;
}

double Body::E_k(void)
{
	double E_t = 0.5*mass*(r_dot.transpose()*r_dot)(0);
	Eigen::Vector3d omega_p = 2.0*caams::L(p)*p_dot;
	double E_r = 0.5*(omega_p.transpose()*J_p*omega_p)(0);
	return E_t + E_r;
}

DatumBody::DatumBody(void):
    Body(
		Eigen::Vector3d::Zero(),
		Eigen::Vector4d(1.0, 0.0, 0.0, 0.0),
		Eigen::Vector3d::Zero(),
		Eigen::Vector4d::Zero(),
        1.0,
		Eigen::Matrix3d::Identity())
{

}

void DatumBody::Draw(void)
{
	glm::dmat3x3 Aglm = E2GLM(caams::Ap(p));
    glm::dmat4x4 Abody(Aglm);
	glm::dvec3 rglm = E2GLM(r);
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
		Eigen::Vector3d r,
		Eigen::Vector4d p,
		Eigen::Vector3d r_dot,
		Eigen::Vector4d p_dot,
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
	glm::dmat3x3 Aglm = E2GLM(caams::Ap(p));
	glm::dmat4x4 Abody(Aglm);
	glm::dvec3 rglm = E2GLM(r);
	glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 radj(0.0,0.0,-length/2.0);
    glm::dmat4x4 Aadj = glm::rotate(M_PI/2.0,glm::dvec3(0.0,1.0,0.0))*
                        glm::translate(radj);
    glm::dmat4x4 Amodel = Atrans*Abody*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
	WireCylinder(radius,length,(GLint)5,(GLint)1);
    glPopMatrix();
}


CylinderYaxis::CylinderYaxis(
		Eigen::Vector3d r,
		Eigen::Vector4d p,
		Eigen::Vector3d r_dot,
		Eigen::Vector4d p_dot,
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
	glm::dmat3x3 Aglm = E2GLM(caams::Ap(p));
	glm::dmat4x4 Abody(Aglm);
	glm::dvec3 rglm = E2GLM(r);
	glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 radj(0.0,0.0,-length/2.0);
    glm::dmat4x4 Aadj = glm::rotate(-M_PI/2.0,glm::dvec3(1.0,0.0,0.0))*
                        glm::translate(radj);
    glm::dmat4x4 Amodel = Atrans*Abody*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
	WireCylinder(radius,length,(GLint)5,(GLint)1);
    glPopMatrix();
}

CylinderZaxis::CylinderZaxis(
		Eigen::Vector3d r,
		Eigen::Vector4d p,
		Eigen::Vector3d r_dot,
		Eigen::Vector4d p_dot,
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
	glm::dmat3x3 Aglm = E2GLM(caams::Ap(p));
	glm::dmat4x4 Abody(Aglm);
	glm::dvec3 rglm = E2GLM(r);
	glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 radj(0.0,0.0,-length/2.0);
    glm::dmat4x4 Aadj = glm::translate(radj);
    glm::dmat4x4 Amodel = Atrans*Abody*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
	WireCylinder(radius,length,(GLint)5,(GLint)1);
    glPopMatrix();
}

Cuboid::Cuboid(
		Eigen::Vector3d r,
		Eigen::Vector4d p,
		Eigen::Vector3d r_dot,
		Eigen::Vector4d p_dot,
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
	glm::dmat3x3 Aglm = E2GLM(caams::Ap(p));
	glm::dmat4x4 Abody(Aglm);
	glm::dvec3 rglm = E2GLM(r);
	glm::dmat4x4 Atrans = glm::translate(rglm);
    glm::dvec3 s(dx,dy,dz);
    glm::dmat4x4 Ascale = glm::scale(s);
    glm::dmat4x4 Amodel = Atrans*Abody*Ascale;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Amodel));
	WireCube(1.0);
    glPopMatrix();

}

Sphere::Sphere(
		Eigen::Vector3d r,
		Eigen::Vector4d p,
		Eigen::Vector3d r_dot,
		Eigen::Vector4d p_dot,
		double mass,
		double radius):
	Body(r,p,r_dot,p_dot,mass,caams::J_p_sphere(mass,radius)),
	radius(radius)
{

}

void Sphere::Draw(void)
{
	glm::dmat3x3 Aglm = E2GLM(caams::Ap(p));
	glm::dmat4x4 Abody(Aglm);
	glm::dvec3 rglm = E2GLM(r);
	glm::dmat4x4 Atrans = glm::translate(rglm);
	glm::dmat4x4 Amodel = Atrans*Abody;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glMultMatrixd(glm::value_ptr(Amodel));
	WireSphere(radius,16,8);
	glPopMatrix();

}


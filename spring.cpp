/*
 * zero length spring/damper between a sphere and ground
 *
 * Equations come from CAAMS equation (11.17)
 *
 */

#include "caams.hpp"
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

double cp_phi = 3.0*M_PI/4.0; // connection point phi
double cp_theta = 0.0; // connection point theta
double mass=1.0;
double radius=0.5;

// rotation about the negative z-axis
caams::matrix p_phi(caams::pAA(cp_phi,caams::matrix(3,1,0.0,0.0,-1.0)));
// rotation about the y-axis
caams::matrix p_theta(caams::pAA(cp_theta,caams::matrix(3,1,0.0,1.0,0.0)));
caams::matrix s0(3,1,0.0,radius,0.0);
caams::matrix A_phi(caams::G(p_phi)*~caams::L(p_phi));
caams::matrix A_theta(caams::G(p_theta)*~caams::L(p_theta));
caams::matrix s_p(A_theta*A_phi*s0); // spring connection point in object coordinates
caams::matrix s_s(s_p); // spring base in world coordinates
caams::matrix g(3,1,0.0,-10.0,0.0); // gravitational vector
caams::matrix J_p(caams::J_p_sphere(mass,radius));
caams::matrix q0(7,1,0.0,0.0,0.0,1.0,0.0,0.0,0.0);
caams::matrix q_dot0(7,1,caams::zeros);

double k_spring = 1.0e4;
double omega_n = std::sqrt(k_spring/mass);
double c_damp = 2.0*mass*omega_n;
double f_n = omega_n/2.0/M_PI;
double dt=1.0/f_n/60.0;

bool paused=true;

caams::matrix system_solve(caams::matrix q, caams::matrix q_dot){
    caams::matrix p(q.sub(4,1,4,1));
    caams::matrix p_dot(q_dot.sub(4,1,4,1));
    caams::matrix r(q.sub(3,1,1,1));
    caams::matrix r_dot(q_dot.sub(3,1,1,1));
    caams::matrix s(caams::G(p)*~caams::L(p)*s_p);
    caams::matrix s_c(r + s);
    caams::matrix d_s(s_s - s_c);
    caams::matrix f_spring(3,1);
    caams::matrix f_damp(3,1);
    double mag_d_s = caams::norm(d_s);
    if(mag_d_s == 0.0){
        f_spring = caams::matrix(3,1,caams::zeros);
        f_damp = caams::matrix(3,1,caams::zeros);
    }else{
        f_spring = k_spring*d_s;
        caams::matrix u_d_s((1.0/mag_d_s)*d_s);
        caams::matrix v_c(r_dot + 2.0*caams::G(p_dot)*~caams::L(p)*s_p);
        caams::matrix v_c_r((u_d_s*~u_d_s)*v_c);
        f_damp = -c_damp*v_c_r;
    }
    caams::matrix f_gravity(mass*g);
    caams::matrix f_total(f_spring+f_damp+f_gravity);
    caams::matrix n_total(caams::SS(s)*(f_spring+f_damp));
    // CAAMS equation 11.25
    caams::matrix twoL(2.0*caams::L(p));
    caams::matrix twoL_dot(2.0*caams::L(p_dot));
    caams::matrix H(~twoL_dot*J_p*twoL);

    caams::matrix y(8,1);
    y.sub(f_total,1,1);
    y.sub(2.0*~caams::G(p)*n_total - 2.0*H*p_dot,4,1);
    y.sub(-~p_dot*p_dot,8,1);

    caams::matrix A(8,8,caams::zeros);
    A.sub(mass*caams::matrix(3,3,caams::init_identity),1,1);
    A.sub(~twoL*J_p*twoL,4,4);
    A.sub(~p,8,4);
    A.sub(p,4,8);

    caams::matrix x(A.inverse()*y);

    caams::matrix q_ddot(x.sub(7,1,1,1));

    return q_ddot;
}

void normalize_p(caams::matrix &q){
    caams::matrix p(q.sub(4,1,4,1));
    p = (1.0/caams::norm(p))*p;
    q.sub(p,4,1);
}

void system_advance(caams::matrix &q,caams::matrix &q_dot, double dt){
    caams::matrix k_q_dot(7,4);
    caams::matrix k_q_ddot(7,4);
    caams::matrix q_norm(7,1);

    k_q_dot.sub(q_dot,1,1);
    k_q_ddot.sub(system_solve(q,k_q_dot.sub(7,1,1,1)),1,1);

    q_norm = q + (dt/2.0)*k_q_dot.sub(7,1,1,1);
    normalize_p(q_norm);
    k_q_dot.sub(q_dot+(dt/2.0)*k_q_ddot.sub(7,1,1,1),1,2);
    k_q_ddot.sub(system_solve(q_norm,k_q_dot.sub(7,1,1,2)),1,2);

    q_norm = q + (dt/2.0)*k_q_dot.sub(7,1,1,2);
    normalize_p(q_norm);
    k_q_dot.sub(q_dot+(dt/2.0)*k_q_ddot.sub(7,1,1,2),1,3);
    k_q_ddot.sub(system_solve(q_norm,k_q_dot.sub(7,1,1,3)),1,3);

    q_norm = q + dt*k_q_dot.sub(7,1,1,3);
    normalize_p(q_norm);
    k_q_dot.sub(q_dot+dt*k_q_ddot.sub(7,1,1,3),1,4);
    k_q_ddot.sub(system_solve(q_norm,k_q_dot.sub(7,1,1,4)),1,4);

    q = q + (dt/6.0)*k_q_dot*caams::matrix(4,1,1.0,2.0,2.0,1.0);
    normalize_p(q);
    q_dot = q_dot + (dt/6.0)*k_q_ddot*caams::matrix(4,1,1.0,2.0,2.0,1.0);
    caams::matrix p(q.sub(4,1,4,1));
    caams::matrix p_dot(q_dot.sub(4,1,4,1));
    caams::matrix sigma(~p_dot*p);
    p_dot -= *sigma.data *p;
    q_dot.sub(p_dot,4,1);
		
}

void render_sphere(caams::matrix &A, caams::matrix &r, double radius){
    glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
    glm::dmat4x4 Amodel(Aglm);
    glm::dvec3 rglm = glm::make_vec3(r.data);
    glm::dmat4x4 Atrans = glm::translate(glm::dmat4x4(1.0),rglm);
    glm::dmat4x4 Aadj = glm::rotate(glm::dmat4(1.0),M_PI/2,glm::dvec3(-1.0,0.0,0.0));
    glm::dmat4x4 Aeff = Atrans*Amodel*Aadj;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(glm::value_ptr(Aeff));
    glutSolidSphere(radius,8,8);
    glPopMatrix();
}

void render_cylinder_x_axis(caams::matrix &A,caams::matrix &r,double radius,double length){
	glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
	glm::dmat4x4 Amodel(Aglm);
	glm::dvec3 rglm = glm::make_vec3(r.data);
	glm::dmat4x4 Atrans = glm::translate(glm::dmat4x4(1.0),rglm);
	glm::dvec3 radj(0.0,0.0,-length/2.0);
	glm::dmat4x4 Aadj = glm::translate(
		glm::rotate(glm::dmat4(1.0),M_PI/2.0,glm::dvec3(0.0,1.0,0.0)),radj);
	Amodel = Atrans*Amodel*Aadj;
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glMultMatrixd(glm::value_ptr(Amodel));
	glutWireCylinder(radius,length,(GLint)5,(GLint)1);
	glPopMatrix();
}

void render_cylinder_y_axis(caams::matrix &A,caams::matrix &r,double radius,double length){
	glm::dmat3x3 Aglm = glm::make_mat3x3((~A).data);
	glm::dmat4x4 Amodel(Aglm);
	glm::dvec3 rglm = glm::make_vec3(r.data);
	glm::dmat4x4 Atrans = glm::translate(glm::dmat4x4(1.0),rglm);
	glm::dvec3 radj(0.0,0.0,-length/2.0);
	glm::dmat4x4 Aadj = glm::translate(
		glm::rotate(glm::dmat4(1.0),-M_PI/2.0,glm::dvec3(1.0,0.0,0.0)),radj);
	Amodel = Atrans*Amodel*Aadj;
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glMultMatrixd(glm::value_ptr(Amodel));
	glutWireCylinder(radius,length,(GLint)5,(GLint)1);
	glPopMatrix();
}

void system_render(void){
    caams::matrix p(q0.sub(4,1,4,1));
    caams::matrix r(q0.sub(3,1,1,1));
    caams::matrix A(caams::G(p)*~caams::L(p));
    r.print("r");

    caams::matrix p_dot(q_dot0.sub(4,1,4,1));
    caams::matrix c(~p_dot*p);
    c.print("velocity contstraint");



    render_sphere(A,r,radius);
}

glm::dmat4 camera_rotation(1.0);
glm::dmat4 camera_translation(glm::translate(glm::dmat4(1.0),glm::dvec3(0.0,0.0,3.0)));

#define NEAR_PLANE 0.1
#define FAR_PLANE 30.0
#define FOV (M_PI*45.0/180.0)

void display(void){
	int window_width;
	int window_height;
	
	window_width = glutGet(GLUT_WINDOW_WIDTH);
	window_height = glutGet(GLUT_WINDOW_HEIGHT);

	glm::dmat4x4 m_projection = glm::perspective(
		FOV,(double)window_width/window_height,NEAR_PLANE,FAR_PLANE);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(glm::value_ptr(m_projection));
	
	glm::dmat4 camera_to_world = camera_translation*camera_rotation;
	glm::dmat4 world_to_camera = glm::inverse(camera_to_world);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(glm::value_ptr(world_to_camera));
	
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	system_render();
	
	glutSwapBuffers();
	
	if(!paused)
        system_advance(q0,q_dot0,dt);
}

#define TIMER_INTERVAL (1000.0/60.0)
double timer_interval=0.0;


void timerFunc( int value )
{
	glutPostRedisplay();
	timer_interval += TIMER_INTERVAL;
	int ti_int = (int)floor(timer_interval);
	glutTimerFunc( ti_int, timerFunc, 0 );
	timer_interval -= ti_int;
}

void keyboardFunc(unsigned char key, int x, int y)
{
	switch(key){
		case ' ':
			paused=paused?false:true;
			break;
		case 27:
			glutLeaveMainLoop();
			break;
		default:
			break;
	}
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowSize (1280, 720);
	glutInitWindowPosition (100, 100);
    glutCreateWindow ("Spring");
	glutDisplayFunc(display);
	glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
	glutKeyboardFunc(keyboardFunc);
	glClearColor(0.0,0.0,0.0,0.0);
	glDrawBuffer(GL_BACK);
	glColor3d(1.0,1.0,1.0);
	glutMainLoop();
	return 0;
}

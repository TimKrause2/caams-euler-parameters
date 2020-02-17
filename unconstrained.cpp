/*
 * Unconstrained rotational motion of a free body
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


caams::matrix J_p_cylinder_x_axis(double m, double r, double l){
	r*=r;
	l*=l;
	double Jxx = m*r/2.0;
	double Jyy = m/12.0*(3*r+l);
	double Jzz = Jyy;
	double i[3][3]=
	{
		{Jxx,0.0,0.0},
		{0.0,Jyy,0.0},
		{0.0,0.0,Jzz}
	};
	return caams::matrix(3,3,(double*)i);
}

caams::matrix J_p_cylinder_y_axis(double m, double r, double h){
	r*=r;
	h*=h;
	double Jxx = m/12.0*(3*r+h);
	double Jyy = m*r/2.0;
	double Jzz = Jxx;
	double i[3][3]=
	{
		{Jxx,0.0,0.0},
		{0.0,Jyy,0.0},
		{0.0,0.0,Jzz}
	};
	return caams::matrix(3,3,(double*)i);
}

double cyl_radius = 0.5;
double cyl_height = 0.05;
double cyl_mass = 1.0;

caams::matrix omega0(3,1, 0.0, 2.0*M_PI, 0.0);
caams::matrix p0(caams::pAA(30.0*M_PI/180.0,caams::matrix(3,1, 0.0,0.0,1.0)));
caams::matrix p_dot0(0.5*~caams::G(p0)*omega0);
caams::matrix Jp(J_p_cylinder_y_axis(cyl_mass, cyl_radius, cyl_height));

double dt=1.0/60.0;
double dt_max=0.001;

bool paused=true;

caams::matrix system_solve(caams::matrix p, caams::matrix p_dot){
    caams::matrix H(4.0*~caams::L(p_dot)*Jp*caams::L(p));
    caams::matrix R(4,4);
    R.sub(caams::L(p)*H,1,1);
    R.sub(~p_dot,4,1);

    caams::matrix y(-R*p_dot);

    caams::matrix A(4,4);
    A.sub(2.0*Jp*caams::L(p),1,1);
    A.sub(~p,4,1);

    caams::matrix x(A.inverse()*y);
	
    return x;
}

void normalize_p(caams::matrix &p){
    p = (1.0/caams::norm(p))*p;
}

double delta_t_recommended(caams::matrix p,caams::matrix p_dot,caams::matrix p_ddot){
    caams::matrix twoLp(2.0*caams::L(p));
    caams::matrix omega_p(twoLp*p_dot);
    caams::matrix alpha_p(twoLp*p_ddot);
	double theta_min=0.1*2.0*M_PI;
    return std::min(theta_min/caams::norm(omega_p),
                    std::sqrt(theta_min/caams::norm(alpha_p)));
}
	
void system_advance(caams::matrix &p,caams::matrix &p_dot, double dt){
    caams::matrix k_p_dot(4,4);
    caams::matrix k_p_ddot(4,4);
    caams::matrix p_norm(4,1);
	double t=0.0;
	while( t<dt ){
		k_p_dot.sub(p_dot,1,1);
        k_p_ddot.sub(system_solve(p,k_p_dot.sub(4,1,1,1)),1,1);
		
        double dt_sim = delta_t_recommended(p,k_p_dot.sub(4,1,1,1),k_p_ddot.sub(4,1,1,1));
		if( dt_sim > dt_max )
			dt_sim = dt_max;
		if( (t+dt_sim) > dt )
			dt_sim = dt-t;
		
        p_norm = p + (dt_sim/2.0)*k_p_dot.sub(4,1,1,1);
		normalize_p(p_norm);
        k_p_dot.sub(p_dot+(dt_sim/2.0)*k_p_ddot.sub(4,1,1,1),1,2);
        k_p_ddot.sub(system_solve(p_norm,k_p_dot.sub(4,1,1,2)),1,2);
		
        p_norm = p + (dt_sim/2.0)*k_p_dot.sub(4,1,1,2);
		normalize_p(p_norm);
        k_p_dot.sub(p_dot+(dt_sim/2.0)*k_p_ddot.sub(4,1,1,2),1,3);
        k_p_ddot.sub(system_solve(p_norm,k_p_dot.sub(4,1,1,3)),1,3);
		
        p_norm = p + dt_sim*k_p_dot.sub(4,1,1,3);
		normalize_p(p_norm);
        k_p_dot.sub(p_dot+dt_sim*k_p_ddot.sub(4,1,1,3),1,4);
        k_p_ddot.sub(system_solve(p_norm,k_p_dot.sub(4,1,1,4)),1,4);
		
		p = p + (dt_sim/6.0)*k_p_dot*caams::matrix(4,1,1.0,2.0,2.0,1.0);
		normalize_p(p);
		p_dot = p_dot + (dt_sim/6.0)*k_p_ddot*caams::matrix(4,1,1.0,2.0,2.0,1.0);
		
		t += dt_sim;
	}
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
    caams::matrix A(caams::G(p0)*~caams::L(p0));
    caams::matrix r(3,1, 0.0, 0.0, 0.0);
    render_cylinder_y_axis(A,r,cyl_radius,cyl_height);
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

    // calculate the angular momentum
    caams::matrix omega(2.0*caams::G(p0)*p_dot0);
    caams::matrix A(caams::G(p0)*~caams::L(p0));
    caams::matrix J(A*Jp*~A);
    caams::matrix L(J*omega);
    L.print("angular momentum");
    caams::matrix E(0.5*~omega*L);
    E.print("rotational energy");

	system_render();
	
	glutSwapBuffers();
	
	if(!paused)
		system_advance(p0,p_dot0,dt);
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
    glutCreateWindow ("Unconstrained");
	glutDisplayFunc(display);
	glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
	glutKeyboardFunc(keyboardFunc);
	glClearColor(0.0,0.0,0.0,0.0);
	glDrawBuffer(GL_BACK);
	glColor3d(1.0,1.0,1.0);
	glutMainLoop();
	return 0;
}

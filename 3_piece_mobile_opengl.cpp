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

caams::matrix sj_jacobian(caams::matrix p, caams::matrix s_p ){
	caams::matrix B(3,7);
	B.sub( caams::matrix(3,3,caams::init_identity), 1, 1);
	B.sub( 2.0*(caams::G(p)*caams::a_minus(s_p)+s_p*~p), 1, 4);
	return B;
}

caams::matrix sj_gamma(caams::matrix p_dot, caams::matrix s_p){
	return -2.0*(caams::G(p_dot)*caams::a_minus(s_p)*p_dot+s_p*(~p_dot*p_dot));
}

bool gyro_flag=false;
double m1=1.0;
double radius1=0.01;
double length1=0.25;
double m2=0.5;
double radius2=0.005;
double length2=0.1;
double radius2_gyro=0.05;
double height2=0.01;
double m3=0.5;
double radius3=0.005;
double length3=0.1;
caams::matrix s_p1_1(3,1, 0.0,radius1*2.0,0.0);
caams::matrix s_p2_1(3,1);
caams::matrix s_p2_2(3,1);
caams::matrix s_p3_1(3,1, length1/2.0,-radius1*2.0,0.0);
caams::matrix s_p3_3(3,1, 0.0,radius3*2.0,0.0);
caams::matrix p1(4,1, 1.0,0.0,0.0,0.0);
caams::matrix p2(caams::pAA(M_PI/4.0,caams::matrix(3,1, 0.0,0.0,1.0)));
//caams::matrix p2(4,1, 1.0,0.0,0.0,0.0);
caams::matrix p3(4,1, 1.0,0.0,0.0,0.0);
caams::matrix omega_p1(3,1, 0.0,0.0,0.0);
//caams::matrix omega_p2(3,1, 0.0,0.0,0.0);
caams::matrix omega_p3(3,1, 0.0,0.0,0.0);
caams::matrix p_dot1(0.5*~caams::L(p1)*omega_p1);
caams::matrix p_dot2(4,1);
caams::matrix p_dot3(0.5*~caams::L(p3)*omega_p3);
caams::matrix g(3,1, 0.0,-9.81,0.0);
caams::matrix f1(m1*g);
caams::matrix f2(m2*g);
caams::matrix f3(m3*g);
double k_omega = -1.0e-5;
caams::matrix k_omega1(3,3);
caams::matrix k_omega2(3,3);
caams::matrix k_omega3(3,3);
caams::matrix J_p1(J_p_cylinder_x_axis(m1,radius1,length1));
caams::matrix J_p2(3,3);
caams::matrix J_p3(J_p_cylinder_x_axis(m3,radius3,length3));

caams::matrix p0(12,1);
caams::matrix p_dot0(12,1);

double dt=1.0/60.0;
double dt_max=0.001;

bool paused=true;

void mobileInit(void){
//	std::cout << "gyro_flag:" << gyro_flag << std::endl;
	
	k_omega1 = caams::matrix(3,3, k_omega*10,0.0,0.0, 0.0,k_omega*40.0,0.0, 0.0,0.0,k_omega*40.0);
	k_omega3 = caams::matrix(3,3, k_omega/4.0,0.0,0.0, 0.0,k_omega,0.0, 0.0,0.0,k_omega);
	
	caams::matrix omega_p2(3,1);
	
	if(gyro_flag){
		s_p2_1 = caams::matrix(3,1, -length1/2.0,radius1,0.0);
		s_p2_2 = caams::matrix(3,1, 0.0,-0.05,0.0);
		J_p2 = J_p_cylinder_y_axis(m2,radius2_gyro,height2);
		k_omega2 = caams::matrix(3,3, k_omega,0.0,0.0, 0.0,k_omega/4.0,0.0, 0.0,0.0,k_omega);
		omega_p2 = caams::matrix(3,1, 0.0,10.0*2.0*M_PI,0.0);
	}else{
		s_p2_1 = caams::matrix(3,1, -length1/2.0,-radius1*2.0,0.0);
		s_p2_2 = caams::matrix(3,1, 0.0,radius2*2.0,0.0);
		J_p2 = J_p_cylinder_x_axis(m2,radius2,length2);
		k_omega2 = caams::matrix(3,3, k_omega/4.0,0.0,0.0, 0.0,k_omega,0.0, 0.0,0.0,k_omega);
		omega_p2 = caams::matrix(3,1, 0.0,0.0,0.0);
	}
	p_dot2 = 0.5*~caams::L(p2)*omega_p2;
}


enum int_mode{
	NEWTON,
	RK4
};

int_mode im=NEWTON;

caams::matrix system_solve(caams::matrix p, caams::matrix p_dot,
						   caams::matrix const &f1, caams::matrix const &f2,
						   caams::matrix const &f3){
	caams::matrix p1(p.sub(4,1,1,1));
	caams::matrix p2(p.sub(4,1,5,1));
	caams::matrix p3(p.sub(4,1,9,1));
	caams::matrix p_dot1(p_dot.sub(4,1,1,1));
	caams::matrix p_dot2(p_dot.sub(4,1,5,1));
	caams::matrix p_dot3(p_dot.sub(4,1,9,1));
	caams::matrix twoLp1(2.0*caams::L(p1));
	caams::matrix twoLp2(2.0*caams::L(p2));
	caams::matrix twoLp3(2.0*caams::L(p3));
	caams::matrix twoLp_dot1(2.0*caams::L(p_dot1));
	caams::matrix twoLp_dot2(2.0*caams::L(p_dot2));
	caams::matrix twoLp_dot3(2.0*caams::L(p_dot3));
	caams::matrix omega_p1(twoLp1*p_dot1);
	caams::matrix omega_p2(twoLp2*p_dot2);
	caams::matrix omega_p3(twoLp3*p_dot3);
	caams::matrix n_p1(k_omega1*omega_p1);
	caams::matrix n_p2(k_omega2*omega_p2);
	caams::matrix n_p3(k_omega3*omega_p3);
	caams::matrix n_star1(~twoLp1*n_p1);
	caams::matrix n_star2(~twoLp2*n_p2);
	caams::matrix n_star3(~twoLp3*n_p3);
	caams::matrix H1(~twoLp_dot1*J_p1*twoLp1);
	caams::matrix H2(~twoLp_dot2*J_p2*twoLp2);
	caams::matrix H3(~twoLp_dot3*J_p3*twoLp3);
	caams::matrix b_star1(2.0*H1*p_dot1);
	caams::matrix b_star2(2.0*H2*p_dot2);
	caams::matrix b_star3(2.0*H3*p_dot3);
	caams::matrix y(33,1);
	y.sub(f1,1,1);
	y.sub(n_star1-b_star1,4,1);
	y.sub(f2,8,1);
	y.sub(n_star2-b_star2,11,1);
	y.sub(f3,15,1);
	y.sub(n_star3-b_star3,18,1);
	y.sub(-~p_dot1*p_dot1,22,1);
	y.sub(-~p_dot2*p_dot2,23,1);
	y.sub(-~p_dot3*p_dot3,24,1);
	y.sub(sj_gamma(p_dot1,s_p1_1),25,1);
	y.sub(sj_gamma(p_dot2,s_p2_2)-sj_gamma(p_dot1,s_p2_1),28,1);
	y.sub(sj_gamma(p_dot3,s_p3_3)-sj_gamma(p_dot1,s_p3_1),31,1);
	
	caams::matrix N1(m1*caams::matrix(3,3,caams::init_identity));
	caams::matrix N2(m2*caams::matrix(3,3,caams::init_identity));
	caams::matrix N3(m3*caams::matrix(3,3,caams::init_identity));
	caams::matrix J_star1(~twoLp1*J_p1*twoLp1);
	caams::matrix J_star2(~twoLp2*J_p2*twoLp2);
	caams::matrix J_star3(~twoLp3*J_p3*twoLp3);
	caams::matrix B11(sj_jacobian(p1,s_p1_1));
	caams::matrix B21(-sj_jacobian(p1,s_p2_1));
	caams::matrix B22(sj_jacobian(p2,s_p2_2));
	caams::matrix B31(-sj_jacobian(p1,s_p3_1));
	caams::matrix B33(sj_jacobian(p3,s_p3_3));
	
	caams::matrix A(33,33,caams::zeros);
	A.sub(N1,1,1);
	A.sub(J_star1,4,4);
	A.sub(N2,8,8);
	A.sub(J_star2,11,11);
	A.sub(N3,15,15);
	A.sub(J_star3,18,18);
	A.sub(p1,4,22);
	A.sub(~p1,22,4);
	A.sub(p2,11,23);
	A.sub(~p2,23,11);
	A.sub(p3,18,24);
	A.sub(~p3,24,18);
	A.sub(B11,25,1);
	A.sub(~B11,1,25);
	A.sub(B21,28,1);
	A.sub(B22,28,8);
	A.sub(~B21,1,28);
	A.sub(~B22,8,28);
	A.sub(B31,31,1);
	A.sub(B33,31,15);
	A.sub(~B31,1,31);
	A.sub(~B33,15,31);
	caams::matrix x(A.inverse()*y);
	
	caams::matrix p_ddot(12,1);
	p_ddot.sub(x.sub(4,1,4,1),1,1);
	p_ddot.sub(x.sub(4,1,11,1),5,1);
	p_ddot.sub(x.sub(4,1,18,1),9,1);
	
	//caams::matrix lambda1(x.sub(3,1,25,1));
	//caams::matrix lambda2(x.sub(3,1,28,1));
	//caams::matrix lambda3(x.sub(3,1,31,1));
	//lambda1.print("lambda1");
	//lambda2.print("lambda2");
	//lambda3.print("lambda3");

	return p_ddot;
}

void normalize_p(caams::matrix &p){
	caams::matrix p1(p.sub(4,1,1,1));
	caams::matrix p2(p.sub(4,1,5,1));
	caams::matrix p3(p.sub(4,1,9,1));
	p.sub((1.0/caams::norm(p1))*p1,1,1);
	p.sub((1.0/caams::norm(p2))*p2,5,1);
	p.sub((1.0/caams::norm(p3))*p3,9,1);
}

double delta_t_recommended(caams::matrix p,caams::matrix p_dot,caams::matrix p_ddot){
	caams::matrix p1(p.sub(4,1,1,1));
	caams::matrix p2(p.sub(4,1,5,1));
	caams::matrix p3(p.sub(4,1,9,1));
	caams::matrix p_dot1(p_dot.sub(4,1,1,1));
	caams::matrix p_dot2(p_dot.sub(4,1,5,1));
	caams::matrix p_dot3(p_dot.sub(4,1,9,1));
	caams::matrix p_ddot1(p_ddot.sub(4,1,1,1));
	caams::matrix p_ddot2(p_ddot.sub(4,1,5,1));
	caams::matrix p_ddot3(p_ddot.sub(4,1,9,1));
	caams::matrix twoLp1(2.0*caams::L(p1));
	caams::matrix twoLp2(2.0*caams::L(p2));
	caams::matrix twoLp3(2.0*caams::L(p3));
	caams::matrix omega_p1(twoLp1*p_dot1);
	caams::matrix omega_p2(twoLp2*p_dot2);
	caams::matrix omega_p3(twoLp3*p_dot3);
	caams::matrix alpha_p1(twoLp1*p_ddot1);
	caams::matrix alpha_p2(twoLp2*p_ddot2);
	caams::matrix alpha_p3(twoLp3*p_ddot3);
	double theta_min=0.1*2.0*M_PI;
	return std::min(
		std::min(theta_min/caams::norm(omega_p1),
		std::min(theta_min/caams::norm(omega_p2),theta_min/caams::norm(omega_p3))),
		std::min(std::sqrt(theta_min/caams::norm(alpha_p1)),
		std::min(std::sqrt(theta_min/caams::norm(alpha_p2)),
				 std::sqrt(theta_min/caams::norm(alpha_p3)))));
}
	
void system_advance(caams::matrix &p,caams::matrix &p_dot, double dt){
	caams::matrix k_p_dot(12,4);
	caams::matrix k_p_ddot(12,4);
	caams::matrix p_norm(12,1);
	double t=0.0;
	while( t<dt ){
		k_p_dot.sub(p_dot,1,1);
		k_p_ddot.sub(system_solve(p,k_p_dot.sub(12,1,1,1),f1,f2,f3),1,1);
		
		double dt_sim = delta_t_recommended(p,k_p_dot.sub(12,1,1,1),k_p_ddot.sub(12,1,1,1));
		if( dt_sim > dt_max )
			dt_sim = dt_max;
		if( (t+dt_sim) > dt )
			dt_sim = dt-t;
		
		p_norm = p + (dt_sim/2.0)*k_p_dot.sub(12,1,1,1);
		normalize_p(p_norm);
		k_p_dot.sub(p_dot+(dt_sim/2.0)*k_p_ddot.sub(12,1,1,1),1,2);
		k_p_ddot.sub(system_solve(p_norm,k_p_dot.sub(12,1,1,2),f1,f2,f3),1,2);
		
		p_norm = p + (dt_sim/2.0)*k_p_dot.sub(12,1,1,2);
		normalize_p(p_norm);
		k_p_dot.sub(p_dot+(dt_sim/2.0)*k_p_ddot.sub(12,1,1,2),1,3);
		k_p_ddot.sub(system_solve(p_norm,k_p_dot.sub(12,1,1,3),f1,f2,f3),1,3);
		
		p_norm = p + dt_sim*k_p_dot.sub(12,1,1,3);
		normalize_p(p_norm);
		k_p_dot.sub(p_dot+dt_sim*k_p_ddot.sub(12,1,1,3),1,4);
		k_p_ddot.sub(system_solve(p_norm,k_p_dot.sub(12,1,1,4),f1,f2,f3),1,4);
		
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
	caams::matrix p1(p0.sub(4,1,1,1));
	caams::matrix p2(p0.sub(4,1,5,1));
	caams::matrix p3(p0.sub(4,1,9,1));
	caams::matrix A1(caams::G(p1)*~caams::L(p1));
	caams::matrix A2(caams::G(p2)*~caams::L(p2));
	caams::matrix A3(caams::G(p3)*~caams::L(p3));
	caams::matrix r1(A1*-s_p1_1);
	caams::matrix r2(r1+A1*s_p2_1-A2*s_p2_2);
	caams::matrix r3(r1+A1*s_p3_1-A3*s_p3_3);
	render_cylinder_x_axis(A1,r1,radius1,length1);
	if(gyro_flag){
		render_cylinder_y_axis(A2,r2,radius2_gyro,height2);
	}else{
		render_cylinder_x_axis(A2,r2,radius2,length2);
	}
	render_cylinder_x_axis(A3,r3,radius3,length3);
}

glm::dmat4 camera_rotation(1.0);
glm::dmat4 camera_translation(glm::translate(glm::dmat4(1.0),glm::dvec3(0.0,0.0,0.5)));

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

char optstring[]="gk:";

struct option longopts[]=
{
	{"gyro",no_argument,NULL,'g'},
	{"k_omega",required_argument,NULL,'k'},
	{NULL,0,NULL,0}
};

int main(int argc, char **argv){
	glutInit(&argc,argv);
	int val;
	int r;
	while( (val=getopt_long(argc,argv,optstring,longopts,NULL)) != -1 ){
		switch(val){
			case 'g':
				gyro_flag=true;
				break;
			case 'k':
				r=sscanf(optarg,"%lg",&k_omega);
				if(r!=1){
					std::cout << "Invalid argument for k_omega." << std::endl;
					std::exit(1);
				}
				break;
			default:
				std::cout << "Invalid arguments." << std::endl;
				std::exit(1);
				break;
		}
	}
	
	mobileInit();
	p0.sub(p1,1,1);
	p0.sub(p2,5,1);
	p0.sub(p3,9,1);
	p_dot0.sub(p_dot1,1,1);
	p_dot0.sub(p_dot2,5,1);
	p_dot0.sub(p_dot3,9,1);
	
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowSize (1280, 720);
	glutInitWindowPosition (100, 100);
	glutCreateWindow ("3 Piece Mobile");
	glutDisplayFunc(display);
	glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
	glutKeyboardFunc(keyboardFunc);
	glClearColor(0.0,0.0,0.0,0.0);
	glDrawBuffer(GL_BACK);
	glColor3d(1.0,1.0,1.0);
	glutMainLoop();
	return 0;
}

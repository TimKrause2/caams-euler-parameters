#include "mbsystem.h"
#include "forces.h"
#include <math.h>
#include <iostream>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>


#define GYRO_RADIUS 0.5
#define GYRO_LENGTH 0.1
#define GYRO_MASS   1.0
#define GYRO_OFFSET 0.75
#define GYRO_THETA  30.0
#define GYRO_OMEGA  10.0*2.0*M_PI

Body *datumBody = new DatumBody;
Body *gyroBody;
Constraint *sphericalJoint;
ForceElement *gravity;

System mbsystem;

void init_system(void)
{
	Eigen::Vector3d s1_p(0.0, -GYRO_OFFSET, 0.0);
	Eigen::Vector3d s2_p(0.0, 0.0, 0.0);

	Eigen::Vector4d p_gyro(caams::pAA(glm::radians(GYRO_THETA),
							 Eigen::Vector3d(0.0, 0.0, 1.0)));
	Eigen::Matrix3d A_gyro(caams::Ap(p_gyro));
	Eigen::Vector3d r_error(0.1, 0.1, 0.0);
	Eigen::Vector3d r_gyro(s2_p - A_gyro*s1_p + r_error);
	Eigen::Vector3d omega_p(0.0, GYRO_OMEGA, 0.0);
	Eigen::Vector3d r_dot_gyro(0.0, 0.0, 0.0);
	Eigen::Vector4d p_dot_gyro(0.5*caams::L(p_gyro).transpose()*omega_p);

    gyroBody = new CylinderYaxis(
                r_gyro,
                p_gyro,
                r_dot_gyro,
                p_dot_gyro,
                GYRO_MASS,
                GYRO_RADIUS,
                GYRO_LENGTH);

    sphericalJoint = new SphericalJoint(gyroBody,datumBody,s1_p,s2_p);

	gravity = new SystemGravityForce(
				Eigen::Vector3d(0.0, -9.81, 0.0),
				mbsystem);

    mbsystem.AddBody(gyroBody);
	mbsystem.AddConstraint(sphericalJoint);
	mbsystem.AddForce(gravity);
    mbsystem.InitializeSolver();
}

glm::dmat4 camera_rotation(1.0);
glm::dmat4 camera_translation(glm::translate(glm::dvec3(0.0,0.0,3.0)));

#define NEAR_PLANE 0.1
#define FAR_PLANE 30.0
#define FOV (M_PI*90.0/180.0)

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

    glColor3f(1.0f,1.0f,1.0f);
    mbsystem.Draw();

    datumBody->Draw();

    glutSwapBuffers();

	mbsystem.Integrate(1/60.0);
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
    init_system();
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
    glutInitWindowSize (1280, 720);
    glutInitWindowPosition (100, 100);
    glutCreateWindow ("mbtest");
    glutDisplayFunc(display);
    glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
    glutKeyboardFunc(keyboardFunc);
    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);
    glColor3d(1.0,1.0,1.0);
    glutMainLoop();
    return 0;
}

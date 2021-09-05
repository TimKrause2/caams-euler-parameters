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
#define GYRO_MASS   0.5
#define GYRO_OFFSET 0.75
#define GYRO_THETA  10.0
#define GYRO_OMEGA  10.0*2.0*M_PI

#define MOBILE1_LENGTH  3.0
#define MOBILE1_RADIUS  0.25
#define MOBILE1_MASS    0.1
#define MOBILE1_THETA   -10.0

#define MOBILE2_LENGTH  1.0
#define MOBILE2_RADIUS  0.1
#define MOBILE2_MASS    0.5
#define MOBILE2_THETA   45.0

Body *datumBody = new DatumBody;
Body *gyroBody;
Body *mobile1Body;
Body *mobile2Body;

Constraint *mobileGyroJoint;
Constraint *mobileDatumJoint;
Constraint *mobileMobileJoint;

ForceElement *gravity;

System mbsystem;

void init_system(void)
{
    caams::matrix s1_p_mobileGyro(3,1,-MOBILE1_LENGTH/2.0,MOBILE1_RADIUS,0.0);
    caams::matrix s2_p_mobileGyro(3,1,0.0,-GYRO_OFFSET,0.0);
    caams::matrix s1_p_mobileDatum(3,1,0.0,0.0,0.0);
    caams::matrix s2_p_mobileDatum(3,1,0.0,MOBILE1_RADIUS+0.5,0.0);
    caams::matrix s1_p_mobileMobile(3,1,MOBILE1_LENGTH/2.0,-MOBILE1_RADIUS,0.0);
    caams::matrix s2_p_mobileMobile(3,1,0.0,MOBILE2_RADIUS+0.15,0.0);

    caams::matrix p_gyro(caams::pAA(glm::radians(GYRO_THETA),
                             caams::matrix(3,1,0.0,0.0,1.0)));
    caams::matrix A_gyro(caams::Ap(p_gyro));
    caams::matrix p_mobile1(caams::pAA(glm::radians(MOBILE1_THETA),
                                       caams::matrix(3,1,0.0,0.0,1.0)));
    caams::matrix A_mobile1(caams::Ap(p_mobile1));
    caams::matrix p_mobile2(caams::pAA(glm::radians(MOBILE2_THETA),
                                       caams::matrix(3,1,0.0,0.0,1.0)));
    caams::matrix A_mobile2(caams::Ap(p_mobile2));

    caams::matrix r_mobile1(s1_p_mobileDatum - A_mobile1*s2_p_mobileDatum);
    caams::matrix r_mobile2(r_mobile1 + A_mobile1*s1_p_mobileMobile - A_mobile2*s2_p_mobileMobile);
    caams::matrix r_gyro(r_mobile1 + A_mobile1*s1_p_mobileGyro - A_gyro*s2_p_mobileGyro);

    caams::matrix omega_p(3,1,0.0,GYRO_OMEGA,0.0);
    caams::matrix r_dot_gyro(3,1,0.0,0.0,0.0);
    caams::matrix p_dot_gyro(0.5*~caams::L(p_gyro)*omega_p);
    caams::matrix r_dot_mobile1(3,1,0.0,0.0,0.0);
    caams::matrix p_dot_mobile1(4,1,0.0,0.0,0.0,0.0);
    caams::matrix r_dot_mobile2(3,1,0.0,0.0,0.0);
    caams::matrix p_dot_mobile2(4,1,0.0,0.0,0.0,0.0);

    gyroBody = new CylinderYaxis(
                r_gyro,
                p_gyro,
                r_dot_gyro,
                p_dot_gyro,
                GYRO_MASS,
                GYRO_RADIUS,
                GYRO_LENGTH);

    mobile1Body = new CylinderXaxis(r_mobile1,p_mobile1,r_dot_mobile1,p_dot_mobile1,
                                    MOBILE1_MASS,MOBILE1_RADIUS,MOBILE1_LENGTH);

    mobile2Body = new CylinderXaxis(r_mobile2,p_mobile2,r_dot_mobile2,p_dot_mobile2,
                                    MOBILE2_MASS,MOBILE2_RADIUS,MOBILE2_LENGTH);

    mobileGyroJoint = new SphericalJoint(mobile1Body,gyroBody,s1_p_mobileGyro,s2_p_mobileGyro);

    mobileDatumJoint = new SphericalJoint(datumBody,mobile1Body,s1_p_mobileDatum,s2_p_mobileDatum);

    mobileMobileJoint = new SphericalJoint(mobile1Body,mobile2Body,s1_p_mobileMobile,s2_p_mobileMobile);

    gravity = new SystemGravityForce(
                caams::matrix(3,1,0.0,-9.81,0.0),
                mbsystem);

    mbsystem.AddBody(gyroBody);
    mbsystem.AddBody(mobile1Body);
    mbsystem.AddBody(mobile2Body);

    mbsystem.AddConstraint(mobileGyroJoint);
    mbsystem.AddConstraint(mobileDatumJoint);
    mbsystem.AddConstraint(mobileMobileJoint);

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

    mbsystem.Integrate(1/(60.0*5.0));
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

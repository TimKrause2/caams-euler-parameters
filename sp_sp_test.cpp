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


#define LENGTH 2.0
#define THETA  5.0
#define BODY_MASS   1.0
#define BODY_RADIUS 0.25
#define BODY_LENGTH 0.5

Body *datumBody = new DatumBody;
Body *pendulumBody;
Constraint *ssJoint;

ForceElement *gravity;

System mbsystem;

void init_system(void)
{
    caams::matrix s1_p(3,1,0.0,-0.5,0.0);
    caams::matrix s2_p(3,1,0.0,0.5,0.0);

    caams::matrix p_body(pAA(glm::radians(THETA),
                             caams::matrix(3,1,0.0,0.0,1.0)));
    caams::matrix A_body(caams::Ap(p_body));
    caams::matrix d(3,1,0.0,-LENGTH,0.0);
    caams::matrix r_body(s2_p + d - A_body*s1_p);
    caams::matrix r_dot_body(3,1,0.0,0.0,0.0);
    caams::matrix p_dot_body(4,1,0.0,0.0,0.0,0.0);

    pendulumBody = new CylinderYaxis(
                r_body,
                p_body,
                r_dot_body,
                p_dot_body,
                BODY_MASS,
                BODY_RADIUS,
                BODY_LENGTH);

    ssJoint = new SphericalSphericalJoint(pendulumBody,datumBody,s1_p,s2_p,LENGTH);

    gravity = new SystemGravityForce(
                caams::matrix(3,1,0.0,-9.81,0.0),
                mbsystem);

    mbsystem.AddBody(pendulumBody);
    mbsystem.AddConstraint(ssJoint);
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

    glutSwapBuffers();

    mbsystem.Integrate(1/600.0);
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
    glutCreateWindow ("SphericalSphericalJoint Test");
    glutDisplayFunc(display);
    glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
    glutKeyboardFunc(keyboardFunc);
    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);
    glColor3d(1.0,1.0,1.0);
    glutMainLoop();
    return 0;
}

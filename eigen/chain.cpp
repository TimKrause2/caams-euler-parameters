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

#define N_LINKS  32

#define DRIVER_MASS       1.0
#define DRIVER_RADIUS     0.5
#define DRIVER_LENGTH     0.1
#define DRIVER_OMEGA      (2.0*M_PI)
#define DRIVER_TORQUE     20.0
#define DRIVER_K_R        (DRIVER_TORQUE/DRIVER_OMEGA)
#define DRIVER_K_T        0.0

#define LINK_MASS         0.2
#define LINK_RADIUS       0.01
#define LINK_LENGTH       0.1
#define LINK_K_T          0.001
#define LINK_K_R          0.003


Body *datumBody = new DatumBody;
Body *driverBody;
Body *linkBodies[N_LINKS];

Constraint *driverSpJoint;
Constraint *driverP1_2Joint;
Constraint *linkSpJoints[N_LINKS];

ForceElement *gravity;
ForceElement *driverFriction;
ForceElement *driverMotor;
ForceElement *linkFriction[N_LINKS];

bool driver_on = false;

System mbsystem;

void init_system(void)
{
	Eigen::Vector3d s1_p_driver(0.0,1.5,0.0);
	Eigen::Vector3d s2_p_driver(0.0,0.0,0.0);
	Eigen::Vector3d s1_p_driverp(0.0,1.0,0.0);
	Eigen::Vector3d s2_p_driverp(0.0,1.0,0.0);
	Eigen::Vector3d s1_p_driver_link(DRIVER_RADIUS,-DRIVER_LENGTH/2.0,0.0);
	Eigen::Vector3d s2_p_link(0.0,LINK_LENGTH/2.0,0.0);
	Eigen::Vector3d s1_p_link_link(0.0,-LINK_LENGTH/2.0,0.0);

	Eigen::Vector3d r_driver(s1_p_driver - s2_p_driver);
	Eigen::Vector4d p(1.0,0.0,0.0,0.0);
	Eigen::Vector4d p_dot = Eigen::Vector4d::Zero();
	Eigen::Vector3d r_dot = Eigen::Vector3d::Zero();

    driverBody = new CylinderYaxis(r_driver,p,r_dot,p_dot,
                                   DRIVER_MASS,DRIVER_RADIUS,DRIVER_LENGTH);

	Eigen::Vector3d r_link(r_driver + s1_p_driver_link - s2_p_link);

    for(int i=0;i<N_LINKS;i++){
        linkBodies[i] = new CylinderYaxis(r_link,p,r_dot,p_dot,
                                          LINK_MASS,LINK_RADIUS,LINK_LENGTH);
        r_link += s1_p_link_link - s2_p_link;
    }

    driverSpJoint = new SphericalJoint(datumBody,driverBody,s1_p_driver,s2_p_driver);
    driverP1_2Joint = new Parallel1_2(datumBody,driverBody,s1_p_driverp,s2_p_driverp);

    linkSpJoints[0] = new SphericalJoint(driverBody,linkBodies[0],s1_p_driver_link,s2_p_link);
    for(int i=1;i<N_LINKS;i++){
        linkSpJoints[i] = new SphericalJoint(linkBodies[i-1],linkBodies[i],s1_p_link_link,s2_p_link);
    }

    gravity = new SystemGravityForce(
				Eigen::Vector3d(0.0,-9.81,0.0),
                mbsystem);
    driverFriction = new BodyDamping(driverBody,DRIVER_K_T,DRIVER_K_R);
    driverMotor = new BodyLocalTorque(driverBody,
									  Eigen::Vector3d::Zero());
    for(int i=0;i<N_LINKS;i++){
        linkFriction[i] = new BodyDamping(linkBodies[i],LINK_K_T,LINK_K_R);
    }

    mbsystem.AddBody(driverBody);
    for(int i=0;i<N_LINKS;i++){
        mbsystem.AddBody(linkBodies[i]);
    }

    mbsystem.AddConstraint(driverSpJoint);
    mbsystem.AddConstraint(driverP1_2Joint);
    for(int i=0;i<N_LINKS;i++){
        mbsystem.AddConstraint(linkSpJoints[i]);
    }

    mbsystem.AddForce(gravity);
    mbsystem.AddForce(driverMotor);
    mbsystem.AddForce(driverFriction);
    for(int i=0;i<N_LINKS;i++){
        mbsystem.AddForce(linkFriction[i]);
    }

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

    BodyLocalTorque *dm = dynamic_cast<BodyLocalTorque*>(driverMotor);
    if(driver_on){
		dm->n = Eigen::Vector3d(0.0,DRIVER_TORQUE,0.0);
    }else{
		dm->n = Eigen::Vector3d::Zero();
    }

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
            if(driver_on)driver_on=false;
            else driver_on=true;
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
	glutCreateWindow ("Chain Demostration");
    glutDisplayFunc(display);
    glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
    glutKeyboardFunc(keyboardFunc);
    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);
    glColor3d(1.0,1.0,1.0);
    glutMainLoop();
    return 0;
}

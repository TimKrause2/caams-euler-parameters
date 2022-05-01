#include "mbsystem.h"
#include "forces.h"
#include "SDL.h"
#include "SDL_joystick.h"
#include "SDL_gamecontroller.h"
#include <math.h>
#include <iostream>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>


#define VEHICLE_MASS         1.0
#define VEHICLE_RADIUS       1.0
#define OMEGA_MAX            (0.5*2.0*M_PI)
#define TORQUE_MAX           10.0
#define K_R                  TORQUE_MAX/OMEGA_MAX
#define VELOCITY_MAX_LONG    50.0
#define VELOCITY_MAX_RADIAL  10.0
#define FORCE_MAX            20.0
#define K_T_LONG             FORCE_MAX/VELOCITY_MAX_LONG
#define K_T_RADIAL           FORCE_MAX/VELOCITY_MAX_RADIAL

#define Z_NEAR 0.1f
#define Z_FAR 1000.0f
#define FOV (M_PI*90.0f/180.0f)

#define ARRAY_DIM 4
#define ARRAY_SIZE 30.0
#define ARRAY_RADIUS 1.0
#define ARRAY_MASS 1.0

Body *vehicleBody;
Body *sceneBodies[ARRAY_DIM][ARRAY_DIM][ARRAY_DIM];
BodyLocalForce *localForce;
BodyLocalTorque *localTorque;
ForceElement *friction;

System mbsystem;

void init_system(void)
{
	Eigen::Vector3d r_vehicle = Eigen::Vector3d::Zero();
	Eigen::Vector4d p_vehicle(1.0, 0.0, 0.0, 0.0);
	Eigen::Vector3d r_dot_vehicle = Eigen::Vector3d::Zero();
	Eigen::Vector4d p_dot_vehicle = Eigen::Vector4d::Zero();

	vehicleBody = new Sphere(
				r_vehicle,
				p_vehicle,
				r_dot_vehicle,
				p_dot_vehicle,
				VEHICLE_MASS,
				VEHICLE_RADIUS);

	localForce = new BodyLocalForce(vehicleBody,Eigen::Vector3d::Zero());
	localTorque = new BodyLocalTorque(vehicleBody,Eigen::Vector3d::Zero());
    Eigen::Matrix3d k_t;
    Eigen::Matrix3d k_r;
    k_t << K_T_RADIAL, 0.0, 0.0,
            0.0, K_T_RADIAL, 0.0,
            0.0, 0.0, K_T_LONG;
    k_r << K_R, 0.0, 0.0,
            0.0, K_R, 0.0,
            0.0, 0.0, K_R;
    friction = new BodyDamping(vehicleBody, k_t, k_r);

	for(int i=0;i<ARRAY_DIM;i++){
		for(int j=0;j<ARRAY_DIM;j++){
			for(int k=0;k<ARRAY_DIM;k++){
				Eigen::Vector3d r_array(
							-ARRAY_SIZE/2.0 + ARRAY_SIZE/(ARRAY_DIM-1)*i,
							-ARRAY_SIZE/2.0 + ARRAY_SIZE/(ARRAY_DIM-1)*j,
							-ARRAY_SIZE/2.0 + ARRAY_SIZE/(ARRAY_DIM-1)*k);
				Eigen::Vector4d p_array(1.0, 0.0, 0.0, 0.0);
				Eigen::Vector3d r_dot_array(0.0, 0.0, 0.0);
				Eigen::Vector4d p_dot_array(0.0, 0.0, 0.0, 0.0);

				sceneBodies[i][j][k] = new Sphere(
							r_array,
							p_array,
							r_dot_array,
							p_dot_array,
							ARRAY_MASS,
							ARRAY_RADIUS);

				mbsystem.AddBody(sceneBodies[i][j][k]);

			}
		}
	}

	mbsystem.AddBody(vehicleBody);
	mbsystem.AddForce(localForce);
	mbsystem.AddForce(localTorque);
	mbsystem.AddForce(friction);
	mbsystem.InitializeSolver();
}

void condition(double &x)
{
	if(x>=0.0){
		x*=x;
	}else{
		x*=-x;
	}
}



int main(int argc, char **argv){
    init_system();

	if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_GAMECONTROLLER)!=0){
		SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
		return 1;
	}

	printf("Platform:%s\n",SDL_GetPlatform());

	printf("Current video driver:%s\n",SDL_GetCurrentVideoDriver());

	SDL_GameController *gc = NULL;
	for(int i=0;i<SDL_NumJoysticks();++i){
		if(SDL_IsGameController(i)){
			gc = SDL_GameControllerOpen(i);
			if(gc){
				break;
			}else{
				std::cout << "Could not open gamecontroller " << i << ": " <<
							 SDL_GetError() << std::endl;
			}
		}
	}

	if(!gc){
		std::cout << "No game controllers found. Exiting." << std::endl;
		return 1;
	}


	SDL_Window *window;
	window = SDL_CreateWindow(
		"Vehicle",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		1280,
		720,
		SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);

	if(window==NULL){
		printf("Couldn't create window: %s\n",SDL_GetError());
		SDL_Quit();
		return 1;
	}

	SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 8 );
	SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 8 );
	SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 8 );
	SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 8 );
	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
	SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24 );
	SDL_GL_SetAttribute( SDL_GL_CONTEXT_PROFILE_MASK,
						 SDL_GL_CONTEXT_PROFILE_COMPATIBILITY );
	SDL_GL_SetAttribute( SDL_GL_CONTEXT_MAJOR_VERSION, 3 );
	SDL_GL_SetAttribute( SDL_GL_CONTEXT_MINOR_VERSION, 1 );

	SDL_GLContext glcontext = SDL_GL_CreateContext(window);

	if(glcontext==NULL){
		printf("Couldn't create GL context: %s\n", SDL_GetError());
		SDL_DestroyWindow(window);
		SDL_Quit();
		return 1;
	}

	if(SDL_GL_SetSwapInterval(1)==-1){
		printf("Swap Interval not supported.\n");
	}

	SDL_ShowWindow(window);

	int width;
	int height;
	SDL_GetWindowSize(window, &width, &height);

	int running=1;
	int gc_zero=0;
	SDL_Event event;
	struct timespec t_last;
	struct timespec t_now;

	short left_x0;
	short left_y0;
	short right_x0;
	short right_y0;
	short left_t0;
	short right_t0;

	clock_gettime( CLOCK_MONOTONIC, &t_last );

	while(running){
		while(SDL_PollEvent(&event)){
			switch(event.type){
				case SDL_KEYDOWN:
					if(event.key.keysym.scancode == SDL_SCANCODE_ESCAPE) running=0;
					break;
				case SDL_QUIT:
					running=0;
					break;
				case SDL_WINDOWEVENT:
					switch(event.window.event){
						case SDL_WINDOWEVENT_RESIZED:
						case SDL_WINDOWEVENT_SIZE_CHANGED:
							width = event.window.data1;
							height = event.window.data2;
							break;
						default:
							break;
					}
				default:
					break;
			}
		}
		glViewport(0, 0, width, height );
		float aspect = (float)width/height;
		glm::mat4x4 Mproj;
		float fov = FOV;
		if(aspect>=1.0){
			Mproj = glm::perspective(fov,aspect,Z_NEAR,Z_FAR);
		}else{
			fov = glm::atan(glm::tan(fov/2.0f)/aspect)*2.0f;
			Mproj = glm::perspective(fov,aspect,Z_NEAR,Z_FAR);
		}
		glMatrixMode(GL_PROJECTION);
		glLoadMatrixf(glm::value_ptr(Mproj));

		glm::dmat3x3 A_rot3 = E2GLM(caams::Ap(vehicleBody->p));
		glm::dmat4x4 A_rot(A_rot3);
		glm::dvec3 rglm = E2GLM(vehicleBody->r);
		glm::dmat4x4 A_trans = glm::translate(rglm);
		glm::dmat4x4 A_cam2wld = A_trans*A_rot;
		glm::dmat4x4 A_wld2cam = glm::inverse(A_cam2wld);

		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixd(glm::value_ptr(A_wld2cam));

		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

		glColor3f(1.0f,1.0f,1.0f);
		mbsystem.Draw();

		SDL_GL_SwapWindow(window);

		clock_gettime( CLOCK_MONOTONIC, &t_now );
		long dt_sec = t_now.tv_sec - t_last.tv_sec;
		long dt_nsec = t_now.tv_nsec - t_last.tv_nsec;
		if(dt_nsec < 0 ){
			dt_nsec += 1000000000;
			dt_sec--;
		}
		double x = dt_sec + (double)dt_nsec*1e-9;

		short left_x = SDL_GameControllerGetAxis(gc, SDL_CONTROLLER_AXIS_LEFTX);
		short left_y = SDL_GameControllerGetAxis(gc, SDL_CONTROLLER_AXIS_LEFTY);
		short right_x = SDL_GameControllerGetAxis(gc, SDL_CONTROLLER_AXIS_RIGHTX);
		short right_y = SDL_GameControllerGetAxis(gc, SDL_CONTROLLER_AXIS_RIGHTY);
		short left_t = SDL_GameControllerGetAxis(gc, SDL_CONTROLLER_AXIS_TRIGGERLEFT);
		short right_t = SDL_GameControllerGetAxis(gc, SDL_CONTROLLER_AXIS_TRIGGERRIGHT);

		std::cout << "Left:<" << left_x << "," << left_y << "> Right:<" <<
					 right_x << "," << right_y << "> Ltrigger:" <<
					 left_t << " Rtrigger:" << right_t << std::endl;
		if(!gc_zero){
			left_x0 = left_x;
			left_y0 = left_y;
			right_x0 = right_x;
			right_y0 = right_y;
			left_t0 = left_t;
			right_t0 = right_t;
			gc_zero = 1;
		}else{
			double lx = (double)(left_x - left_x0)/32768.0;
			double ly = (double)(left_y - left_y0)/32768.0;
			double rx = (double)(right_x - right_x0)/32768.0;
			double ry = (double)(right_y - right_y0)/32768.0;
			double lt = (double)(left_t - left_t0)/32768.0;
			double rt = (double)(right_t - right_t0)/32768.0;

			condition(lx);
			condition(ly);
			condition(rx);
			condition(ry);
			condition(lt);
			condition(rt);

			localForce->f = Eigen::Vector3d(lx*FORCE_MAX,
											0.0,
											ly*FORCE_MAX);
			localTorque->n = Eigen::Vector3d(ry*TORQUE_MAX,
											 -rx*TORQUE_MAX,
											 (lt-rt)*TORQUE_MAX);

		}


		mbsystem.Integrate(x);

		t_last = t_now;
	}
	SDL_GL_DeleteContext(glcontext);
	SDL_DestroyWindow(window);
	SDL_Quit();

    return 0;
}

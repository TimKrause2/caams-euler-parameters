#include "caams.hpp"
#include <math.h>

namespace caams
{
	Eigen::Matrix3x4d G(Eigen::Vector4d const & p){
		Eigen::Matrix3x4d r;
		r << -p(1), p(0),-p(3), p(2),
			 -p(2), p(3), p(0),-p(1),
			 -p(3),-p(2), p(1), p(0);
		return r;
	}
	
	Eigen::Matrix3x4d L(Eigen::Vector4d const & p){
		Eigen::Matrix3x4d r;
		r << -p(1), p(0), p(3),-p(2),
			 -p(2),-p(3), p(0), p(1),
			 -p(3), p(2),-p(1), p(0);
		return r;
	}
	
	Eigen::Matrix3d SS(Eigen::Vector3d const & a){
		Eigen::Matrix3d r;
		r << 0.0 ,-a(2), a(1),
			 a(2), 0.0 ,-a(0),
			-a(1), a(0), 0.0;
		return r;
	}
	
	Eigen::Vector4d pAA(double angle, Eigen::Vector3d const & axis){
		Eigen::Vector4d r;
		angle/=2.0;
		r(0) = cos(angle);
		r.segment<3>(1) = axis.normalized()*sin(angle);
		return r;
	}

	Eigen::Matrix3d Ap(Eigen::Vector4d const & p){
		return G(p)*L(p).transpose();
    }
	
	Eigen::Matrix4d a_plus(Eigen::Vector3d const & a){
		Eigen::Matrix4d r;
		r <<
			 0.0 ,-a(0),-a(1),-a(2),
			 a(0), 0.0 ,-a(2), a(1),
			 a(1), a(2), 0.0 ,-a(0),
			 a(2),-a(1), a(0), 0.0;
		return r;
	}

	Eigen::Matrix4d a_minus(Eigen::Vector3d const & a){
		Eigen::Matrix4d r;
		r <<
			 0.0 ,-a(0),-a(1),-a(2),
			 a(0), 0.0 , a(2),-a(1),
			 a(1),-a(2), 0.0 , a(0),
			 a(2), a(1),-a(0), 0.0;
		return r;
	}

	Eigen::Matrix3d J_p_cylinder_x_axis(double m, double r, double l){
        r*=r;
        l*=l;
        double Jxx = m*r/2.0;
        double Jyy = m/12.0*(3*r+l);
        double Jzz = Jyy;
		Eigen::Matrix3d res;
		res<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return res;
    }

	Eigen::Matrix3d J_p_cylinder_y_axis(double m, double r, double h){
        r*=r;
        h*=h;
        double Jxx = m/12.0*(3*r+h);
        double Jyy = m*r/2.0;
        double Jzz = Jxx;
		Eigen::Matrix3d res;
		res<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return res;
	}

	Eigen::Matrix3d J_p_cylinder_z_axis(double m, double r, double h){
        r*=r;
        h*=h;
        double Jxx = m/12.0*(3*r+h);
        double Jzz = m*r/2.0;
        double Jyy = Jxx;
		Eigen::Matrix3d res;
		res<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return res;
	}

	Eigen::Matrix3d J_p_sphere(double mass, double radius){
        double I_sphere = 2.0/5.0*mass*radius*radius;
		Eigen::Matrix3d r;
		r<<
			I_sphere, 0.0, 0.0,
			0.0, I_sphere, 0.0,
			0.0, 0.0, I_sphere;
		return r;
	}

	Eigen::Matrix3d J_p_cuboid(double mass, double dx, double dy, double dz){
        dx*=dx;
        dy*=dy;
        dz*=dz;
        mass*=1.0/12.0;
        double Jxx = mass*(dy+dz);
        double Jyy = mass*(dx+dz);
        double Jzz = mass*(dx+dy);
		Eigen::Matrix3d r;
		r<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return r;
	}
}

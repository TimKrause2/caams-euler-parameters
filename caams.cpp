#include "caams.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>

namespace caams
{
	matrix G(matrix const & p){
		if(p.rows!=4 || p.cols!=1){
			std::cout << "G : matrix is wrong dimensions" << std::endl;
            assert(0);
		}
		double *s = p.data;
		double i[3][4]=
		{
			{-s[1],s[0],-s[3],s[2]},
			{-s[2],s[3],s[0],-s[1]},
			{-s[3],-s[2],s[1],s[0]}
		};
		matrix r(3,4,(double *)i);
		return r;
	}
	
	matrix L(matrix const & p){
		if(p.rows!=4 || p.cols!=1){
			std::cout << "L : matrix is wrong dimensions" << std::endl;
            assert(0);
		}
		double *s = p.data;
		double i[3][4]=
		{
			{-s[1],s[0],s[3],-s[2]},
			{-s[2],-s[3],s[0],s[1]},
			{-s[3],s[2],-s[1],s[0]}
		};
		matrix r(3,4,(double *)i);
		return r;
	}
	
	matrix SS(matrix const & a){
		if(a.rows!=3 || a.cols!=1){
			std::cout << "SS : matrix is wrong dimensions" << std::endl;
            assert(0);
		}
		double *s = a.data;
		double i[3][3]=
		{
			{0.0,-s[2],s[1]},
			{s[2],0.0,-s[0]},
			{-s[1],s[0],0.0}
		};
		matrix r(3,3,(double *)i);
		return r;
	}
	
	matrix pAA(double angle, matrix const & axis){
		if(!(axis.rows==3 && axis.cols==1)){
			std::cout << "pAA : axis matrix wrong dimensions" << std::endl;
            assert(0);
		}
		matrix r(4,1);
		angle/=2.0;
		r.sub( std::cos(angle), 1, 1);
		r.sub( (std::sin(angle)/norm(axis))*axis, 2, 1);
		return r;
	}

    matrix Ap(matrix const & p){
        if(!(p.rows==4 && p.cols==1)){
            std::cout << "Ap : p wrong dimensions." << std::endl;
            assert(0);
        }
        return G(p)*~L(p);
    }
	
	matrix a_plus(matrix const & a){
		if(!(a.rows==3 && a.cols==1)){
			std::cout << "a_plus : matrix is wrong dimensions" << std::endl;
            assert(0);
		}
		double *d = a.data;
		double i[4][4]=
		{
			{0.0,-d[0],-d[1],-d[2]},
			{d[0],0.0,-d[2],d[1]},
			{d[1],d[2],0.0,-d[0]},
			{d[2],-d[1],d[0],0.0}
		};
		matrix r(4,4,(double*)i);
		return r;
	}
	
	matrix a_minus(matrix const & a){
		if(!(a.rows==3 && a.cols==1)){
			std::cout << "a_minus : matrix is wrong dimensions" << std::endl;
            assert(0);
		}
		double *d = a.data;
		double i[4][4]=
		{
			{0.0,-d[0],-d[1],-d[2]},
			{d[0],0.0,d[2],-d[1]},
			{d[1],-d[2],0.0,d[0]},
			{d[2],d[1],-d[0],0.0}
		};
		matrix r(4,4,(double*)i);
		return r;
	}
	
	double norm(matrix const & a){
		if(!(a.rows==1 || a.cols==1)){
			std::cout << "norm : not a vector" << std::endl;
            assert(0);
		}
		double *s = a.data;
		double r=0.0;
		for(long i=a.rows*a.cols;i;i--){
			r += *s * *(s++);
		}
		return std::sqrt(r);
	}
	
	matrix x_axis(matrix const & a){
		if(a.rows!=3 || a.cols!=3){
			std::cout << "x_axis : not a 3 by 3 matrix" << std::endl;
            assert(0);
		}
		matrix x(3,1);
		double *d = x.data;
		double *s = a.data;
		for(long i=3;i;i--){
			*(d++) = *s;
			s+=a.cols;
		}
		return x;
	}
	
	matrix y_axis(matrix const & a){
		if(a.rows!=3 || a.cols!=3){
			std::cout << "y_axis : not a 3 by 3 matrix" << std::endl;
            assert(0);
		}
		matrix y(3,1);
		double *d = y.data;
		double *s = &a.data[1];
		for(long i=3;i;i--){
			*(d++) = *s;
			s+=a.cols;
		}
		return y;
	}
	
	matrix z_axis(matrix const & a){
		if(a.rows!=3 || a.cols!=3){
			std::cout << "z_axis : not a 3 by 3 matrix" << std::endl;
            assert(0);
		}
		matrix z(3,1);
		double *d = z.data;
		double *s = &a.data[2];
		for(long i=3;i;i--){
			*(d++) = *s;
			s+=a.cols;
		}
		return z;
	}
	
    matrix normalize(matrix const & a){
        if(!(a.rows==1||a.cols==1)){
            std::cout << "normalize : not a vector" << std::endl;
            assert(0);
        }
        return (1.0/norm(a))*a;
    }

    long maxComponent(matrix const &a){
        if(!(a.rows==1||a.cols==1)){
            std::cout << "maxComponent : not a vector." << std::endl;
            assert(0);
        }
        long i_max = 0;
        double d_max = fabs(a.data[0]);
        long N = a.rows*a.cols;
        for(long i_test=1;i_test<N;i_test++){
            double d_test = fabs(a.data[i_test]);
            if(d_test>d_max){
                d_max = d_test;
                i_max = i_test;
            }
        }
        i_max++;
        return i_max;
    }

    matrix J_p_cylinder_x_axis(double m, double r, double l){
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

    matrix J_p_cylinder_y_axis(double m, double r, double h){
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

    matrix J_p_cylinder_z_axis(double m, double r, double h){
        r*=r;
        h*=h;
        double Jxx = m/12.0*(3*r+h);
        double Jzz = m*r/2.0;
        double Jyy = Jxx;
        double i[3][3]=
        {
            {Jxx,0.0,0.0},
            {0.0,Jyy,0.0},
            {0.0,0.0,Jzz}
        };
        return caams::matrix(3,3,(double*)i);
    }

    matrix J_p_sphere(double mass, double radius){
        double I_sphere = 2.0/5.0*mass*radius*radius;
        return caams::matrix(3,3, I_sphere,0.0,0.0, 0.0,I_sphere,0.0, 0.0,0.0,I_sphere);
    }

    matrix J_p_cuboid(double mass, double dx, double dy, double dz){
        dx*=dx;
        dy*=dy;
        dz*=dz;
        mass*=1.0/12.0;
        double Jxx = mass*(dy+dz);
        double Jyy = mass*(dx+dz);
        double Jzz = mass*(dx+dy);
        return matrix(3,3, Jxx,0.0,0.0, 0.0,Jyy,0.0, 0.0,0.0,Jzz);
    }
}

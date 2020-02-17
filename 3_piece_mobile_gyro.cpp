#include "caams.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

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

caams::matrix J_p_cylinder_z_axis(double m, double r, double h){
	r*=r;
	h*=h;
	double Jxx = m/12.0*(3*r+h);
	double Jyy = Jxx;
	double Jzz = m*r/2.0;
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

double m1=1.0;
double radius1=0.01;
double length1=0.25;
double m2=0.5;
double radius2=0.35;
double height2=0.05;
double m3=0.5;
double radius3=0.005;
double length3=0.1;
caams::matrix s_p1_1(3,1, 0.0,0.0,radius1*2.0);
caams::matrix s_p2_1(3,1, -length1/2.0,0.0,-radius1*2.0);
caams::matrix s_p2_2(3,1, 0.0,0.0,-0.3);
caams::matrix s_p3_1(3,1, length1/2.0,0.0,-radius1*2.0);
caams::matrix s_p3_3(3,1, 0.0,0.0,radius3*2.0);
caams::matrix p1(4,1, 1.0,0.0,0.0,0.0);
caams::matrix p2(caams::pAA(M_PI/4.0,caams::matrix(3,1, 0.0,1.0,0.0)));
//caams::matrix p2(4,1, 1.0,0.0,0.0,0.0);
caams::matrix p3(4,1, 1.0,0.0,0.0,0.0);
caams::matrix omega_p1(3,1, 0.0,0.0,0.0);
caams::matrix omega_p2(3,1, 0.0,0.0,20.0*M_PI);
caams::matrix omega_p3(3,1, 0.0,0.0,0.0);
caams::matrix p_dot1(0.5*~L(p1)*omega_p1);
caams::matrix p_dot2(0.5*~L(p2)*omega_p2);
caams::matrix p_dot3(0.5*~L(p3)*omega_p3);
caams::matrix g(3,1, 0.0,0.0,-9.81);
caams::matrix f1(m1*g);
caams::matrix f2(m2*g);
caams::matrix f3(m3*g);
double k_omega = 0.0;
caams::matrix k_omega1(3,3, k_omega,0.0,0.0, 0.0,k_omega*4.0,0.0, 0.0,0.0,k_omega*4.0);
caams::matrix k_omega2(3,3, k_omega/4.0,0.0,0.0, 0.0,k_omega,0.0, 0.0,0.0,k_omega);
caams::matrix k_omega3(k_omega2);
caams::matrix J_p1(J_p_cylinder_x_axis(m1,radius1,length1));
caams::matrix J_p2(J_p_cylinder_z_axis(m2,radius2,height2));
caams::matrix J_p3(J_p_cylinder_x_axis(m3,radius3,length3));

double dt=1.0/30.0;
double dt_max=0.001;

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

void print_angles(caams::matrix &p){
	caams::matrix p1(p.sub(4,1,1,1));
	caams::matrix p2(p.sub(4,1,5,1));
	caams::matrix p3(p.sub(4,1,9,1));
	double theta1 = std::acos(p1.data[0])*2*180.0/M_PI*(std::signbit(p1.data[2])?-1.0:1.0);
	double theta2 = std::acos(p2.data[0])*2*180.0/M_PI*(std::signbit(p2.data[2])?-1.0:1.0);
	double theta3 = std::acos(p3.data[0])*2*180.0/M_PI*(std::signbit(p3.data[2])?-1.0:1.0);
	std::cout << "theta1:" << theta1 << std::endl;
	std::cout << "theta2:" << theta2 << std::endl;
	std::cout << "theta3:" << theta3 << std::endl;
}

void print_z_axis(caams::matrix &p){
	caams::matrix p1(p.sub(4,1,1,1));
	caams::matrix p2(p.sub(4,1,5,1));
	caams::matrix p3(p.sub(4,1,9,1));
	caams::matrix A1(caams::G(p1)*~caams::L(p1));
	caams::matrix A2(caams::G(p2)*~caams::L(p2));
	caams::matrix A3(caams::G(p3)*~caams::L(p3));
	caams::matrix z1(caams::z_axis(A1));
	caams::matrix z2(caams::z_axis(A2));
	caams::matrix z3(caams::z_axis(A3));
	z1.print("z1");
	z2.print("z2");
	z3.print("z3");
	std::cout << std::endl;
}

int main(int narg, char **argc){
	caams::matrix p0(12,1);
	p0.sub(p1,1,1);
	p0.sub(p2,5,1);
	p0.sub(p3,9,1);
	caams::matrix p_dot0(12,1);
	p_dot0.sub(p_dot1,1,1);
	p_dot0.sub(p_dot2,5,1);
	p_dot0.sub(p_dot3,9,1);

	double t=0.0;
	
	while(t<1.0){
		print_z_axis(p0);
		system_advance(p0,p_dot0,dt);
		t+=dt;
	}
	print_z_axis(p0);
	return 0;
}

#include "caams.hpp"
#include <cmath>

caams::matrix J_p_cylinder(double m, double r, double l){
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

caams::matrix sj_jacobian(caams::matrix p, caams::matrix s_p ){
	caams::matrix B(3,7);
	B.sub( caams::matrix(3,3,caams::init_identity), 1, 1);
	B.sub( 2.0*(caams::G(p)*caams::a_minus(s_p)+s_p*~p), 1, 4);
	return B;
}

caams::matrix sj_gamma(caams::matrix p_dot, caams::matrix s_p){
    double *de = p_dot.data;
    double *s = s_p.data;

    return -2.0*caams::matrix(3,1,
        de[0]*(de[2]*s[3-1]-de[3]*s[2-1]+2.0*de[0]*s[1-1])
        + de[1]*(de[3]*s[3-1]+de[2]*s[2-1]+2.0*de[1]*s[1-1])
        + de[2]*(de[0]*s[3-1]+de[1]*s[2-1])
        + de[3]*(de[1]*s[3-1]-de[0]*s[2-1]),
        de[0]*(-de[1]*s[3-1]+2.0*de[0]*s[2-1]+de[3]*s[1-1])
        + de[1]*(de[2]*s[1-1]-de[0]*s[3-1])
        + de[2]*(de[3]*s[3-1]+2.0*de[2]*s[2-1]+de[1]*s[1-1])
        + de[3]*(de[2]*s[3-1]+de[0]*s[1-1]),
        de[0]*(2.0*de[0]*s[3-1]+de[1]*s[2-1]-de[2]*s[1-1])
        + de[1]*(de[0]*s[2-1]+de[3]*s[1-1])
        + de[2]*(de[3]*s[2-1]-de[0]*s[1-1])
        + de[3]*(2.0*de[3]*s[3-1]+de[2]*s[2-1]+de[1]*s[1-1]));
}

double m1=1.0;
double radius1=0.01;
double length1=0.25;
double m2=0.5;
double radius2=0.005;
double length2=0.1;
double m3=0.5;
double radius3=0.005;
double length3=0.1;
caams::matrix s_p1_1(3,1, 0.0,0.0,radius1*2.0);
caams::matrix s_p2_1(3,1, -length1/2.0,0.0,-radius1*2.0);
caams::matrix s_p2_2(3,1, 0.0,0.0,radius2*2.0);
caams::matrix s_p3_1(3,1, length1/2.0,0.0,-radius1*2.0);
caams::matrix s_p3_3(3,1, 0.0,0.0,radius3*2.0);
caams::matrix p1(4,1, 1.0,0.0,0.0,0.0);
//caams::matrix p2(pAA(M_PI/4.0,caams::matrix(3,1, 0.0,1.0,0.0)));
caams::matrix p2(4,1, 1.0,0.0,0.0,0.0);
caams::matrix p3(4,1, 1.0,0.0,0.0,0.0);
caams::matrix omega_p1(3,1, 0.0,0.0,0.0);
caams::matrix omega_p2(3,1, 0.0,0.0,0.0);
caams::matrix omega_p3(3,1, 0.0,0.0,0.0);
caams::matrix p_dot1(0.5*~L(p1)*omega_p1);
caams::matrix p_dot2(0.5*~L(p2)*omega_p2);
caams::matrix p_dot3(0.5*~L(p3)*omega_p3);
caams::matrix g(3,1, 0.0,0.0,-9.81);
caams::matrix f1(m1*g);
caams::matrix f2(m2*g);
caams::matrix f3(m3*g);
double k_omega = -0.1;
caams::matrix k_omega1(3,3, k_omega,0.0,0.0, 0.0,k_omega*4.0,0.0, 0.0,0.0,k_omega*4.0);
caams::matrix k_omega2(3,3, k_omega/4.0,0.0,0.0, 0.0,k_omega,0.0, 0.0,0.0,k_omega);
caams::matrix k_omega3(k_omega2);
caams::matrix J_p1(J_p_cylinder(m1,radius1,length1));
caams::matrix J_p2(J_p_cylinder(m2,radius2,length2));
caams::matrix J_p3(J_p_cylinder(m3,radius3,length3));


caams::matrix system_solve(caams::matrix &p, caams::matrix &p_dot,
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
	
	caams::matrix lambda1(x.sub(3,1,25,1));
	caams::matrix lambda2(x.sub(3,1,28,1));
	caams::matrix lambda3(x.sub(3,1,31,1));
	lambda1.print("lambda1");
	lambda2.print("lambda2");
	lambda3.print("lambda3");

	return p_ddot;
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
	
	caams::matrix p_ddot(system_solve(p0,p_dot0,f1,f2,f3));
	p_ddot.print("p_ddot");
	return 0;
}

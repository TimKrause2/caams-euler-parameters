#include "caams.hpp"
#include "math.h"

int main(int narg, char **argc){
	caams::matrix eye4x4(4,4,caams::init_identity);
	eye4x4.print("eye4x4");
	
	double i[4][4]=
	{
		{1.0,2.0,3.0,4.0},
		{2.0,1.0,4.0,5.0},
		{3.0,2.0,1.0,4.0},
		{4.0,1.0,2.0,3.0}
	};
	caams::matrix m44(4,4,(double*)i);
	m44.print("m44");
	caams::matrix m44_inverse(m44.inverse());
	m44_inverse.print("m44_inverse");
	caams::matrix m44_test( m44*m44_inverse );
	m44_test.print("m44_test");

	caams::matrix axis(3,1,0.0,0.0,10.0);
	
	caams::matrix p(pAA(M_PI/4.0,axis));
	
	p.print("p");
	
	caams::matrix A(G(p)*~L(p));
	A.print("A");






	return 0;
}
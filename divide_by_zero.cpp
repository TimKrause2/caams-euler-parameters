#include <cmath>
#include <iostream>

int main(int narg, char **argv){
	double r=1.0/0.0;
	std::cout << "r:" << r << std::endl;
	std::cout << "min:" << std::min(1.0,std::min(1.0/0.0,std::sqrt(1.0/0.0))) << std::endl;
	return 0;
}

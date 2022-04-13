#ifndef FORCE_ELEMENT_H
#define FORCE_ELEMENT_H

#include "caams.hpp"

class ForceElement
{
public:
    ForceElement(void){}
    ~ForceElement(void){}
	void AccumulateForce(Eigen::VectorXd &y_rhs,
						 Eigen::Vector3d f,
                         long row){
		y_rhs.segment<3>(row) += f;
    }
	void AccumulateTorque(Eigen::VectorXd &y_rhs,
						  Eigen::Vector4d n_star,
                          long row){
		y_rhs.segment<4>(row) += n_star;
	}
	virtual void Apply(Eigen::VectorXd &y_rhs)=0;
    virtual void Draw(void){}
};

#endif

#ifndef FORCE_ELEMENT_H
#define FORCE_ELEMENT_H

#include "caams.hpp"


class ForceElement
{
public:
    ForceElement(void){}
    ~ForceElement(void){}
    void AccumulateForce(caams::matrix &y_rhs,
                         caams::matrix f,
                         long row){
        caams::matrix f0(y_rhs.sub(3,1,row,1));
        f0+=f;
        y_rhs.sub(f0,row,1);

    }
    void AccumulateTorque(caams::matrix &y_rhs,
                          caams::matrix n,
                          long row){
        caams::matrix n0(y_rhs.sub(4,1,row,1));
        n0+=n;
        y_rhs.sub(n0,row,1);
    }
    virtual void Apply(caams::matrix &y_rhs)=0;
    virtual void Draw(void){}
};

#endif

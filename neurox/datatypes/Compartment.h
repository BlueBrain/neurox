#pragma once

#include "neurox/neurox.h"

#include <memory>
#include <vector>

using namespace std;

class Compartment
{
  public:
    Compartment(){};
    ~Compartment(){};

    void setSolverValues(double a, double b, double d, double v, double rhs, double area)
    {
        this->a=a; this->b=b; this->d=d; this->v=v; this->rhs=rhs; this->area=area;
    };

    void addChild(shared_ptr<Compartment> child)
    {
        shared_ptr<Compartment> branch = left==nullptr ? left : right;
        branch=child;
    };

  private:
    double a,b,d,v,rhs,area;

    vector<int> mechanismsIds;
    vector<double> data;
    vector<Datum> pdata;

    shared_ptr<Compartment> left, right;
};

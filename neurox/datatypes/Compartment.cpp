#pragma once

#include "neurox/neurox.h"

#include <memory>
#include <vector>

using namespace std;

Compartment::Compartment(){};
Compartment::~Compartment(){};

void Compartment::setSolverValues(double a, double b, double d, double v, double rhs, double area)
{
    this->a=a; this->b=b; this->d=d; this->v=v; this->rhs=rhs; this->area=area;
};

void Compartment::addChild(shared_ptr<Compartment> child)
{
    assert(left==nullptr || right==nullptr);
    shared_ptr<Compartment> branch = left==nullptr ? left : right;
    branch=child;
};

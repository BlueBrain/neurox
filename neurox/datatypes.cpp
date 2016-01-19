#pragma once

#include "neurox/datatypes.h"

Branch::Branch():
    pdata(NULL), data(NULL), children(NULL), childrenCount(0), 
    n(0), a(NULL), b(NULL), d(NULL), rhs(NULL),
    pdata(NULL), data(NULL), children(NULL)
{
}

Branch::~Branch()
{
    delete [] b; b=NULL;
    delete [] d; d=NULL;
    delete [] a; a=NULL;
    delete [] rhs; rhs=NULL;
    delete [] v; v=NULL;
    delete [] area; area=NULL;

    delete [] pdata; pdata=NULL;
    delete [] children; children=NULL;
}

void Branch::serialize(byte *& bytes_out, int & size_out)
{

}

void Branch::deserialize(byte * bytes_in, int & size_in)
{

}


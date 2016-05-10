#include "neurox/neurox.h"

Branch::Branch():
    pdata(NULL), data(NULL), children(NULL), childrenCount(0), 
    n(0), a(NULL), b(NULL), d(NULL), rhs(NULL)
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
    //TODO missing deletes here
}

void Branch::serialize(byte *& bytes_out, int & size_out)
{

}

void Branch::deserialize(const byte * bytes_in, const int size_in)
{

}


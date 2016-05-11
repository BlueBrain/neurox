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

hpx_action_t Branch::initialize = 0;
int Branch::initialize_handler(const byte * branch_serial_input,  const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t branch_addr = hpx_thread_current_target();
    Branch * branch = NULL;
    if (!hpx_gas_try_pin(branch_addr, (void**) &branch))
        return HPX_RESEND;

    //do the work

    //unpin and return success
    hpx_gas_unpin(branch_addr);
    return HPX_SUCCESS;
}

void Branch::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_POINTER, HPX_SIZE_T);
}

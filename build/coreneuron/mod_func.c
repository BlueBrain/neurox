#include <stdio.h>
extern int nrnmpi_myid;
extern int nrn_nobanner_;
extern int ;

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, " Additional mechanisms from files\n");


    fprintf(stderr, "\n\n");
  }


}

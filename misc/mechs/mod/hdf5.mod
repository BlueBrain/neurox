COMMENT
If the local variable step method is used then the only variables that should
be added are variables of the cell in which this FileRecord
has been instantiated.
ENDCOMMENT

NEURON {
	POINT_PROCESS HDF5
}

PARAMETER {
}

INITIAL {
}


VERBATIM

//#include "/home/fschuerm/hdf5/include/hdf5.h"
//#include "/opt/hdf5/include/hdf5.h"
#include <hdf5.h>

extern double* hoc_pgetarg(int iarg);
extern double* getarg(int iarg);
extern char* gargstr(int iarg);
extern int hoc_is_str_arg(int iarg);
extern int nrnmpi_numprocs;
extern int nrnmpi_myid;
extern int ifarg(int iarg);
extern double chkarg(int iarg, double low, double high);

hid_t oFile, pdataset, sdataset;
herr_t status;
hid_t pdataspace, sdataspace;
int np, pdim;
int ns, sdim;

int pcount, scount;
double* pdata; /* x, y, z, d */
int* sdata; /* pointref, type, parent */

ENDVERBATIM

CONSTRUCTOR {
VERBATIM {


  ns = 0; np = 0;
  pdim = 0;
  sdim = 0;

  pcount = 0; 
  scount = 0; 

}
ENDVERBATIM
}

DESTRUCTOR {
VERBATIM {
  H5Sclose(pdataspace);
  H5Dclose(pdataset);
  H5Sclose(sdataspace);
  H5Dclose(sdataset);
  H5Fclose(oFile);

  free(pdata);
  free(sdata);
}
ENDVERBATIM
}

PROCEDURE setup() { 
VERBATIM { 
  if (hoc_is_str_arg(1)) {
    oFile = H5Fcreate(gargstr(1), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
  } else {
    oFile = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
  }

  if (ifarg(2)) {
     ns = (int) *getarg(2);
  }

  if (ifarg(3)) {
    sdim =  (int) *getarg(3);
  }

  if (ifarg(4)) {
     np = (int) *getarg(4);
  }

  if (ifarg(5)) {
     pdim = (int) *getarg(5);
  }


  pdata = (double *) malloc(np*sizeof(double)*pdim); 
  sdata = (int *) malloc(ns*sizeof(int)*sdim); 

  hsize_t pdims[2]; pdims[0] = np; pdims[1] = pdim;
  hsize_t sdims[2]; sdims[0] = ns; sdims[1] = sdim;
  

  pdataspace = H5Screate_simple(2, pdims, NULL);
  pdataset = H5Dcreate1(oFile, "/points", H5T_NATIVE_DOUBLE, pdataspace, H5P_DEFAULT);
  sdataspace = H5Screate_simple(2, sdims, NULL);
  sdataset = H5Dcreate1(oFile, "/structure", H5T_NATIVE_INT, sdataspace, H5P_DEFAULT);

}
ENDVERBATIM
}

PROCEDURE addSection() { 
VERBATIM { 
  int i;
  for (i=0; i<sdim; i++) {
    sdata[scount*sdim + i]   = (int) *getarg(i+1);
  }
  scount = scount + 1;
}
ENDVERBATIM
}

PROCEDURE addPoint() { 
VERBATIM { 
  int i;
  for (i=0; i<pdim; i++) {
    pdata[pcount*pdim + i]   = *getarg(i+1);
  }
  pcount = pcount + 1;
}
ENDVERBATIM
}

PROCEDURE save() { 
VERBATIM { 

  status = H5Dwrite(pdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pdata);
  status = H5Dwrite(sdataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdata);
	
}
ENDVERBATIM
}

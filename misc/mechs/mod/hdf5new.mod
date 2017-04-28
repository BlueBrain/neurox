COMMENT
ENDCOMMENT

NEURON {
	POINT_PROCESS HDF5NEW
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

hid_t neuronGroup, structureGroup, rawGroup;
hid_t oFile, pdataset, sdataset, tdataset;
herr_t status;
hid_t pdataspace, sdataspace, tdataspace;
int np, pdim;
int ns, sdim;
int nt, tdim;

int pcount, scount, tcount;
double* pdata; /* x, y, z, d */
int* sdata; /* pointref, parent */
int* tdata; /* type - in case of raw morphology has same magnitude as sdata */

ENDVERBATIM

CONSTRUCTOR {
VERBATIM {


  ns = 0; np = 0; nt = 0;
  pdim = 0;
  sdim = 0;
  tdim = 0;

  pcount = 0; 
  scount = 0; 
  tcount = 0;

}
ENDVERBATIM
}

DESTRUCTOR {
VERBATIM {
  H5Sclose(pdataspace);
  H5Dclose(pdataset);
  H5Sclose(sdataspace);
  H5Dclose(sdataset);
  H5Sclose(tdataspace);
  H5Dclose(tdataset);
  H5Fclose(oFile);

  free(pdata);
  free(sdata);
  free(tdata);
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
     nt = ns;
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
  tdim = 1;


  pdata = (double *) malloc(np*sizeof(double)*pdim); 
  sdata = (int *) malloc(ns*sizeof(int)*sdim); 
  tdata = (int *) malloc(nt*sizeof(int)*tdim); 

  hsize_t pdims[2]; pdims[0] = np; pdims[1] = pdim;
  hsize_t sdims[2]; sdims[0] = ns; sdims[1] = sdim;
  hsize_t tdims[1]; tdims[0] = nt; 
  

  neuronGroup = H5Gcreate1(oFile, "/neuron1", 0);
  rawGroup = H5Gcreate1(neuronGroup, "raw", 0);
  structureGroup = H5Gcreate1(neuronGroup, "structure", 0);

  pdataspace = H5Screate_simple(2, pdims, NULL);
  pdataset = H5Dcreate1(rawGroup, "points", H5T_NATIVE_DOUBLE, pdataspace, H5P_DEFAULT);
  sdataspace = H5Screate_simple(2, sdims, NULL);
  sdataset = H5Dcreate1(structureGroup, "raw", H5T_NATIVE_INT, sdataspace, H5P_DEFAULT);
  tdataspace = H5Screate_simple(1, tdims, NULL);
  tdataset = H5Dcreate1(structureGroup, "sectiontype", H5T_NATIVE_INT, tdataspace, H5P_DEFAULT);
  
  // attribute 
  hid_t atype, attr;
  hid_t astring;
  atype  = H5Screate(H5S_SCALAR);
  astring = H5Tcopy(H5T_C_S1); H5Tset_size(astring, 10);
  attr = H5Acreate1(neuronGroup, "creator", astring, atype, H5P_DEFAULT);
  H5Awrite(attr, astring, "neuronHDF5"); 
  H5Aclose(attr); 

  int version = 1;
  atype = H5Screate(H5S_SCALAR);
  attr  = H5Acreate1(neuronGroup, "version", H5T_NATIVE_INT, atype, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, &version);
  H5Aclose(attr); 	



}
ENDVERBATIM
}

PROCEDURE addSection() { 
VERBATIM { 
  sdata[scount*sdim]   = (int) *getarg(1);
  tdata[tcount*tdim]   = (int) *getarg(2);
  sdata[scount*sdim+1]   = (int) *getarg(3);
  scount = scount + 1;
  tcount = tcount + 1;
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
  status = H5Dwrite(tdataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tdata);
	
}
ENDVERBATIM
}

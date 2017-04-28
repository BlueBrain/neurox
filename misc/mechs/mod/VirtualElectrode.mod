COMMENT
/**
 * @file VirtualElectrode.mod
 * @brief 
 * @author reimann
 * @date 2010-12-30
 * @remark Copyright © BBP/EPFL 2005-2011; All rights reserved. Do not distribute without further notice.
 */
ENDCOMMENT

VERBATIM

#include <complex.h>

#include "fftw3.h"

typedef struct {
double *real_data;

fftw_complex *freq_data;

double *time_values;
double *voltage_values;
int num_values;
int N;
int fillPointer;
int replayPointer;
double amplitude;
int isActive;
} signal;

typedef struct {
double phase;
double scale1;
double scale2;
double freq;
} parameters;


ENDVERBATIM

NEURON {
	ARTIFICIAL_CELL VirtualElectrode
	POINTER ptr
	POINTER param_ptr
}
ASSIGNED{
	ptr
	param_ptr
}

CONSTRUCTOR{
	VERBATIM

	signal** tempTable = (signal**)(&(_p_ptr));
	signal* tbl = 0;
	tbl = (signal*)hoc_Emalloc(sizeof(signal));
	tbl->num_values = 0;
	tbl->N = 0;					
	*tempTable = tbl;
	/*Only phase parameter is used actually at this stage. Electrode is modeled as a constant phase element.
	//In the future though I will model the freqency dependent attenuation better, once i have better data*/
	parameters ** tempParam = (parameters**)(&(_p_param_ptr));
	parameters *params = 0;
	params = (parameters*)hoc_Emalloc(sizeof(parameters));
	params->phase = 0.6;
	params->scale1 = 0.004;
	params->scale2 = 0.00001125;
	params->freq = 10000; 
	
	*tempTable = tbl;
	*tempParam = params;

	ENDVERBATIM
}

PROCEDURE setLgth(){
VERBATIM
if(ifarg(1)){
int nVals = (int)*getarg(1);
signal ** tempSig = (signal**)(&(_p_ptr));
signal *sig = (signal*) *tempSig;
sig->num_values = nVals;
sig->time_values = (double*)hoc_Emalloc(sizeof(double)*nVals);
sig->voltage_values = (double*)hoc_Emalloc(sizeof(double)*nVals);
sig->fillPointer = 0;
sig->amplitude = 0;
}
return;
ENDVERBATIM
}

PROCEDURE addPoint(){
VERBATIM
if(ifarg(2)){
	double tVal = (double)*getarg(1);
	double vVal = (double)*getarg(2);
	signal ** tempSig = (signal**)(&(_p_ptr));
	signal *sig = (signal*) *tempSig;
	if(sig->fillPointer < sig->num_values){
		sig->time_values[sig->fillPointer] = tVal;
		sig->voltage_values[sig->fillPointer++] = vVal;
		if(abs(vVal) > sig->amplitude)
			sig->amplitude = vVal;
	} else {
		printf("Already full!\n");
	}
}
return;
ENDVERBATIM
}

PROCEDURE convert(){
VERBATIM
//printf("Starting conversion process!\n");
int i,j;
signal ** tempSig = (signal**)(&(_p_ptr));
signal *sig = (signal*) *tempSig;
parameters ** tempParam = (parameters**)(&(_p_param_ptr));
parameters *params = (parameters*) *tempParam;
int tLgth = (int)(sig->time_values[sig->fillPointer-1])*10;
sig->N = tLgth;
sig->real_data = (double*) hoc_Emalloc(sizeof(double)*sig->N);

sig->freq_data = (fftw_complex*) hoc_Emalloc(sizeof(fftw_complex)*(sig->N/2+1));

double fillVal = 0;
int index = 0;
int tVecPtr = 0;
int multFac = 10;
while(index < sig->N){
	if(sig->time_values[tVecPtr]*multFac <= index){
		fillVal = sig->voltage_values[tVecPtr++];
	}	
	sig->real_data[index++] = fillVal;
}
sig->replayPointer = 0;
double *frequency = (double*)hoc_Emalloc(sizeof(double)*(sig->N/2+1));

double dd = 0;
for(i = 0; i < (sig->N/2+1); ++i){
		frequency[i] = (5000)*dd/(sig->N/2);
		dd+=1;
}

fftw_plan p;
p = fftw_plan_dft_r2c_1d(sig->N, sig->real_data, sig->freq_data, FFTW_ESTIMATE);
fftw_execute(p);
fftw_destroy_plan(p);

for(i = 1; i < (sig->N/2+1); ++i){		
	fftw_complex freqFac = cpow((frequency[i]*I),params->phase);
	/*Here I implement division of complex number myself because.... no reason really*/
	double* foo = (double*) &(sig->freq_data[i]);
	double* bar = (double*) &freqFac;
	foo[0] = (foo[0]*bar[0] + foo[1]*bar[1])/(bar[0]*bar[0]+bar[1]*bar[1]);
	foo[1] = (foo[1]*bar[0]-foo[0]*bar[1])/(bar[0]*bar[0]+bar[1]*bar[1]);
}
/*for(i = (sig->N/2+1); i < sig->N; ++i){	
	sig->freq_data[i] = sig->freq_data[sig->N-i];
}*/
p = fftw_plan_dft_c2r_1d(sig->N, sig->freq_data, sig->real_data, FFTW_ESTIMATE);
fftw_execute(p);
fftw_destroy_plan(p);
for(i = 0; i < (sig->N); ++i){
	sig->real_data[i] = 50*sig->real_data[i]/sig->N;
}
sig->fillPointer = 0;
sig->isActive = 0;

ENDVERBATIM
}


FUNCTION getValue(){
VERBATIM
double distance = 1;
if(ifarg(1)){
	distance = (double)*getarg(1);
}
signal **tempData = (signal**)(&_p_ptr);
signal *sig = (signal*) *tempData;
//printf("%i\n",sig->fillPointer);

	if(sig->fillPointer < sig->num_values && sig->time_values[sig->fillPointer]*10 <= sig->replayPointer){
		sig->fillPointer++;
		if(sig->fillPointer < sig->num_values && sig->voltage_values[sig->fillPointer] != 0){
			sig->isActive = 30;
		}
	}	
	
	
if(sig->replayPointer < sig->N && sig->isActive){
	sig->isActive--;
	return (sig->real_data[sig->replayPointer++])/distance;	
}else{
    sig->replayPointer++;
	return 0;
}
ENDVERBATIM
}

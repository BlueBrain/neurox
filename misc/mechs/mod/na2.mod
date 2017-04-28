:Comment :
:Reference :Traub et. al. J. Neurophysiol, 2003, 89: 909-921

: Adapted by Werner Van Geit @ BBP, 2015 (with help from M.Hines):
: channel detects TTX concentration set by TTXDynamicsSwitch.mod
NEURON {
	SUFFIX na2
	USEION na READ ena WRITE ina
	USEION ttx READ ttxo, ttxi VALENCE 1
	RANGE gna2bar, gna2, ina 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gna2bar = 0.00001 (S/cm2) 
}

ASSIGNED {
	ttxo (mM)
	ttxi (mM)
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gna2	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gna2 = gna2bar*m*m*m*h
	ina = gna2*(v-ena)
}

DERIVATIVE states	{
	if (ttxi == 0.015625 && ttxo > 1e-12) {
		mInf = 0.0
		mTau = 1e-12
		hInf = 1.0
		hTau = 1e-12
	} else {
		rates()
	}
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	if (ttxi == 0.015625 && ttxo > 1e-12) {
		mInf = 0.0
		mTau = 1e-12
		hInf = 1.0
		hTau = 1e-12
	} else {
		rates()
	}
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	UNITSOFF
		mInf = 1/(1+exp((-v-34.5)/10))
        if(v<-26.5){
                mTau = 0.025+0.14*exp((v+26.5)/10)
        }else{
                mTau = 0.02+0.145*exp((-v-26.5)/10)
        }
		hInf = 1/(1+exp((v+59.4)/10.7))
		hTau = 0.15+1.15/(1+exp((v+33.5)/15))
	UNITSON
}

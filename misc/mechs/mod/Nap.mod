:Comment :
:Reference :Magistretti.J and Alonso.A, J. Gen. Physiol, 1999

: Adapted by Werner Van Geit @ BBP, 2015 (with help from M.Hines):
: channel detects TTX concentration set by TTXDynamicsSwitch.mod
NEURON {
	SUFFIX Nap
	USEION na READ ena WRITE ina
	USEION ttx READ ttxo, ttxi VALENCE 1
	RANGE gNapbar, gNap, ina 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNapbar = 0.00001 (S/cm2) 
}

ASSIGNED {
	ttxo (mM)
	ttxi (mM)
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNap	(S/cm2)
	mInf
	mTau
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNap = gNapbar*m*m*m
	ina = gNap*(v-ena)
}

DERIVATIVE states	{
	if (ttxi == 0.015625 && ttxo > 1e-12) {
		mInf = 0.0
		mTau = 1e-12
	} else {
		rates()
	}
	m' = (mInf-m)/mTau
}

INITIAL{
	if (ttxi == 0.015625 && ttxo > 1e-12) {
		mInf = 0.0
		mTau = 1e-12
	} else {
		rates()
	}
	m = mInf
}

PROCEDURE rates(){
	UNITSOFF
		mInf = 1.0000/(1+ exp((v - -44.0000)/-4.8500))
		mTau = 1
	UNITSON
}

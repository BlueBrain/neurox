:Comment :
:Reference :Kinetics were fit to data from Huguenard et al. (1988) and Hamill et al. (1991)

: Adapted by Werner Van Geit @ BBP, 2015 (with help from M.Hines):
: channel detects TTX concentration set by TTXDynamicsSwitch.mod
NEURON {
	SUFFIX Na
	USEION na READ ena WRITE ina
	USEION ttx READ ttxo, ttxi VALENCE 1
	RANGE gNabar, gNa, ina 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNabar = 0.00001 (S/cm2) 
}

ASSIGNED {
	ttxo (mM)
	ttxi (mM)
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNa	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNa = gNabar*m*m*m*h
	ina = gNa*(v-ena)
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
        if(v == -35){
            v = v+0.0001
        }
		mAlpha = (0.182 * (v- -35))/(1-(exp(-(v- -35)/9)))
		mBeta  = (0.124 * (-v -35))/(1-(exp(-(-v -35)/9)))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = 1/(mAlpha + mBeta)
        if(v == -50){
            v = v + 0.0001
        }
		hAlpha = (0.024 * (v- -50))/(1-(exp(-(v- -50)/5)))
        if(v == -75){
            v = v+0.0001
        }
		hBeta  = (0.0091 * (-v - 75))/(1-(exp(-(-v - 75)/5)))
		hInf = 1.0/(1+exp((v- -65)/6.2))
		hTau = 1/(hAlpha + hBeta)
	UNITSON
}

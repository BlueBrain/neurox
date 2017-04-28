TITLE T-calcium channel
: T-type calcium channel for Mala Shah


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER {
	v (mV)
	celsius = 25	(degC)
	gcatbar=.003 (mho/cm2)
	cai = 50.e-6 (mM)
	cao = 2 (mM)
	q10 = 5
	mmin=0.2
	hmin=10
	a0h =0.015
	zetah = 3.5
	vhalfh = -75
	gmh=0.6
	a0m =0.04
	zetam = 2
	vhalfm = -28
	gmm=0.61
	vhm=-60
	vhh=-85
}


NEURON {
	SUFFIX cat_hip
	USEION ca READ cai,cao WRITE ica
        RANGE gcatbar, ica, gcat
        GLOBAL hinf,minf,mtau,htau
}

STATE {
	m h
}

ASSIGNED {
	ica (mA/cm2)
        gcat (mho/cm2)
	hinf
	htau
	minf
	mtau
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcat = gcatbar*m*m*h
	ica = gcat*ghk(v,cai,cao)

}

DERIVATIVE states {	: exact when v held constant
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alph(v(mV)) {
  alph = exp(0.0378*zetah*(v-vhalfh))
}

FUNCTION beth(v(mV)) {
  beth = exp(0.0378*zetah*gmh*(v-vhalfh))
}

FUNCTION alpmt(v(mV)) {
  alpmt = exp(0.0378*zetam*(v-vhalfm))
}

FUNCTION betmt(v(mV)) {
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm))
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a,b, qt
        qt=q10^((celsius-25)/10)
	minf = (1/(1 + exp(-(v-vhm)/10)))
	mtau = betmt(v)/(qt*a0m*(1+alpmt(v)))
	if (mtau<mmin) {mtau=mmin}

	hinf = (1/(1 + exp((v-vhh)/10)))
	htau = beth(v)/(a0h*(1+alph(v)))
	if (htau<hmin) {htau=hmin}
}

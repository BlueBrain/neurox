:Migliore file Modify by Maciej Lazarewicz (mailto:mlazarew@gmu.edu) May/16/2001

TITLE Borg-Graham type generic K-AHP channel

NEURON {
	SUFFIX KahpM95_hip
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,ik, gkahp
        GLOBAL inf,tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	celsius = 6.3	(degC)
	gbar	= .003 	(mho/cm2)
        n	= 4
        cai	= 50.e-6 (mM)
        a0	= 1e8 (/ms-mM-mM-mM-mM)		:b0/(20e-4^4)
        b0	= .5e-2  (/ms)			:0.5/(0.100e3)
        v       	 (mV)
        ek      	 (mV)
	q10=3
}

STATE {	w }

ASSIGNED {
	ik 		(mA/cm2)
        gkahp  		(mho/cm2)
        inf
        tau
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gkahp = gbar*w
	ik = gkahp*(v-ek)
}

INITIAL {
	rate(cai)
	w=inf
}

FUNCTION alp(cai (mM)) {
  alp = a0*cai^n
}

DERIVATIVE state {     : exact when v held constant; integrates over dt step
        rate(cai)
        w' = (inf - w)/tau
}

PROCEDURE rate(cai (mM)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alp(cai)
        tau = 1/(qt*(a + b0))
        inf = a*tau*qt
}

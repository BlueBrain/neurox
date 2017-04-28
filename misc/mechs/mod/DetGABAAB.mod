COMMENT
/**
 * @file DetGABAAB.mod
 * @brief Adapted from ProbGABAA_EMS.mod by Eilif, Michael and Srikanth
 * @author chindemi
 * @date 2014-05-25
 * @remark Copyright © BBP/EPFL 2005-2014; All rights reserved. Do not distribute without further notice.
 */
ENDCOMMENT


TITLE GABAA and GABAB receptor with presynaptic short-term plasticity


COMMENT
GABAA and GABAB receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al. 2002, deterministic
version.
ENDCOMMENT


NEURON {
    THREADSAFE

    POINT_PROCESS DetGABAAB
    RANGE tau_r_GABAA, tau_d_GABAA, tau_r_GABAB, tau_d_GABAB
    RANGE Use, u, Dep, Fac, u0, GABAB_ratio
    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB
    NONSPECIFIC_CURRENT i
    RANGE synapseID, verboseLevel
}


PARAMETER {
    tau_r_GABAA  = 0.2   (ms) : dual-exponential conductance profile
    tau_d_GABAA  = 8     (ms) : IMPORTANT: tau_r < tau_d
    tau_r_GABAB  = 3.5   (ms) : dual-exponential conductance profile :Placeholder value from hippocampal recordings SR
    tau_d_GABAB  = 260.9 (ms) : IMPORTANT: tau_r < tau_d  :Placeholder value from hippocampal recordings
    Use          = 1.0   (1)  : Utilization of synaptic efficacy
    Dep          = 100   (ms) : relaxation time constant from depression
    Fac          = 10    (ms) :  relaxation time constant from facilitation
    e_GABAA      = -80   (mV) : GABAA reversal potential
    e_GABAB      = -97   (mV) : GABAB reversal potential
    gmax         = .001  (uS) : weight conversion factor (from nS to uS)
    u0           = 0          :initial value of u, which is the running value of release probability
    synapseID    = 0
    verboseLevel = 0
    GABAB_ratio  = 0     (1)  : The ratio of GABAB to GABAA
}


ASSIGNED {
    v (mV)
    i (nA)
    i_GABAA (nA)
    i_GABAB (nA)
    g_GABAA (uS)
    g_GABAB (uS)
    g (uS)
    factor_GABAA
    factor_GABAB
}

STATE {
    A_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAA
    B_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAA
    A_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAB
    B_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAB
}


INITIAL{
    LOCAL tp_GABAA, tp_GABAB

    A_GABAA = 0
    B_GABAA = 0

    A_GABAB = 0
    B_GABAB = 0

    tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA) :time to peak of the conductance
    tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) :time to peak of the conductance

    factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) :GABAA Normalization factor - so that when t = tp_GABAA, gsyn = gpeak
    factor_GABAA = 1/factor_GABAA

    factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) :GABAB Normalization factor - so that when t = tp_GABAB, gsyn = gpeak
    factor_GABAB = 1/factor_GABAB
}


BREAKPOINT {
    SOLVE state METHOD cnexp
    g_GABAA = gmax*(B_GABAA-A_GABAA) :compute time varying conductance as the difference of state variables B_GABAA and A_GABAA
    g_GABAB = gmax*(B_GABAB-A_GABAB) :compute time varying conductance as the difference of state variables B_GABAB and A_GABAB
    g = g_GABAA + g_GABAB
    i_GABAA = g_GABAA*(v-e_GABAA) :compute the GABAA driving force based on the time varying conductance, membrane potential, and GABAA reversal
    i_GABAB = g_GABAB*(v-e_GABAB) :compute the GABAB driving force based on the time varying conductance, membrane potential, and GABAB reversal
    i = i_GABAA + i_GABAB
}


DERIVATIVE state{
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
    A_GABAB' = -A_GABAB/tau_r_GABAB
    B_GABAB' = -B_GABAB/tau_d_GABAB
}


NET_RECEIVE (weight,weight_GABAA, weight_GABAB, R, Pr, u, tsyn (ms)){
    LOCAL result
    weight_GABAA = weight
    weight_GABAB = weight * GABAB_ratio

    INITIAL{
            R=1
            u=u0
            tsyn=t
    }

    : calc u at event-
    if (Fac > 0) {
        u = u*exp(-(t - tsyn)/Fac) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    } else {
        u = Use
    }
    if(Fac > 0){
        u = u + Use*(1-u) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    }

    R  = 1 - (1-R) * exp(-(t-tsyn)/Dep) :Probability R for a vesicle to be available for release, analogous to the pool of synaptic
                                        :resources available for release in the deterministic model. Eq. 3 in Fuhrmann et al.
    Pr  = u * R                         :Pr is calculated as R * u (running value of Use)
    R  = R - u * R                      :update R as per Eq. 3 in Fuhrmann et al.

    if( verboseLevel > 0 ) {
        printf("Synapse %f at time %g: R = %g Pr = %g erand = %g\n", synapseID, t, R, Pr, result )
    }

    tsyn = t

    A_GABAA = A_GABAA + Pr*weight_GABAA*factor_GABAA
    B_GABAA = B_GABAA + Pr*weight_GABAA*factor_GABAA
    A_GABAB = A_GABAB + Pr*weight_GABAB*factor_GABAB
    B_GABAB = B_GABAB + Pr*weight_GABAB*factor_GABAB

    if( verboseLevel > 0 ) {
        printf( " vals %g %g %g %g\n", A_GABAA, weight_GABAA, factor_GABAA, weight )
    }
}


FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}

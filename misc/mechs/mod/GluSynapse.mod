COMMENT
/**
 * @file GluSynapse.mod
 * @brief Two state deterministic model of a Glutamatergic Synapse
 * @author chindemi, king
 * @date 2015-05-16
 * @version 0.2.1
 * @acknowledgement Michael Hines, for reviewing many important aspects of model
                    implementation, optimization, and all the technical support. 
                    Francesco Cremonesi, for helping with DE simplification.
 * @remark Copyright © BBP/EPFL 2005-2015; All rights reserved.
           Do not distribute without further notice.
 */
ENDCOMMENT


TITLE Glutamatergic Synapse


COMMENT
In brief:
    Two state deterministic model of a Glutamatergic Synapse.
             ----------------------------------
             |                                |
             |                                v
    Functional Synapse             Potential Synapse
             ^                                |
             |                                |
             ----------------------------------

Functional Synapse components:
    ## AMPAR ##
        Adapted from "ProbAMPANMDA_EMS.mod"
        Suffix: _AMPA

    ## NMDAR ##
        Adapted from "ProbAMPANMDA_EMS.mod"
        Suffix: _NMDA

    ## VDCC ##
        None

    ## Postsynaptic Ca2+ dynamics ##
        See (Graupner and Brunel, 2012)
        Suffix: _GB

    ## Release ##
        Adapted from "ProbAMPANMDA_EMS.mod"
        See (Fuhrmann et al., 2002; Jahr and Stevens, 1990)

    ## LTP/LTD ##
        See (Graupner and Brunel, 2012)
        Suffix: _GB

    ## Elimination ##
        Test
        Suffix: _SE

Potential Synapse components:
    ## Synaptogenesis  ##
        Test
        Suffix: _SG

ENDCOMMENT


NEURON {
    THREADSAFE : TODO Check this keyword
    POINT_PROCESS GluSynapse

    : AMPAR range variables
    RANGE tau_r_AMPA, tau_d_AMPA, i_AMPA, g_AMPA

    : NMDAR range variables
    RANGE tau_r_NMDA, tau_d_NMDA, i_NMDA, g_NMDA
    RANGE mg, mggate

    : Postsynaptic Ca2+ dynamics range variables
    RANGE tau_ca_GB, C_pre_GB, C_post_GB, D_GB

    : Release range variables
    RANGE Use, u, Dep, Fac, u0, Rstate, tsyn, tsyn_fac, Psurv
    POINTER rng_rel

    : LTP/LTD range variables
    RANGE tau_GB, rho_star_GB
    RANGE gamma_p_GB, gamma_d_GB, theta_p_GB, theta_d_GB
    RANGE sigma_GB
    RANGE w0_GB, w1_GB
    RANGE val_rng_GB, Theta_d_GB, Theta_p_GB
    POINTER rng_GB

    : Elimination range variables
    RANGE tau_SE, theta_e_SE

    : Synaptogenesis range variables
    RANGE tau_SG, theta_g_SG

    : Shared range variables
    RANGE g, e, NMDA_ratio, vv, baseweight
    RANGE weight_AMPA, weight_NMDA, factor_AMPA, factor_NMDA

    : Other range variables
    RANGE synapseID, verboseLevel, LTPlasticity, rewiring, synapseState

    : Misc
    NONSPECIFIC_CURRENT i
}


UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (molar) = (1/liter)
    (mM)    = (millimolar)
    (uS)    = (microsiemens)
    (pC)    = (picocoulomb)
}


PARAMETER {
    : AMPAR parameters
    tau_r_AMPA   = 0.2   (ms)   : Dual-exponential conductance profile
    tau_d_AMPA   = 1.7   (ms)   : IMPORTANT: tau_r < tau_d

    : NMDAR paramenters
    tau_r_NMDA   = 0.29  (ms)   : Dual-exponential conductance profile
    tau_d_NMDA   = 43    (ms)   : IMPORTANT: tau_r < tau_d
    alpha_ca_NMDA= 0.2          : TEMPORARY VALUE
    mg           = 1     (mM)   : Initial concentration of mg2+

    : Postsynaptic Ca2+ dynamics parameters
    tau_ca_GB    = 17.63 (ms)   : See SI (Graupner and Brunel, 2012)
    C_pre_GB     = 4.77         : See SI (Graupner and Brunel, 2012)
    C_post_GB    = 1.03         : See SI (Graupner and Brunel, 2012)
    D_GB         = 22.39 (ms)   : See SI (Graupner and Brunel, 2012)

    : Release parameters, just initial values! Use,
    : Dep and Fac are overwritten by BlueBuilder assigned values
    Use          = 1.0   (1)    : Utilization of synaptic efficacy
    Dep          = 100   (ms)   : Relaxation time constant from depression
    Fac          = 10    (ms)   : Relaxation time constant from facilitation
    u0           = 0            : Initial value of u, which is the running value
                                : of release probability

    : LTP/LTD parameters
    w0_GB        = 0.5   (nS)
    w1_GB        = 1.5   (nS)
    rho_star_GB  = 0.5
    tau_GB       = 1413.92  (s)
    theta_d_GB   = 1.0
    theta_p_GB   = 1.3
    gamma_d_GB   = 1321.88
    gamma_p_GB   = 1689.96
    sigma_GB     = 0.82

    : Elimination
    tau_SE       = 1.0   (s)
    theta_e_SE   = 0.2
    factor_SE    = 0.5

    : Synaptogenesis
    tau_SG       = 1.0   (s)
    theta_g_SG   = 0.4
    step_SG      = 0.01

    : Shared parameters
    e            = 0     (mV)   : AMPA and NMDA reversal potential
    gmax         = .001  (uS)   : Weight conversion factor (from nS to uS)
    baseweight   = 0.5   (nS)
    NMDA_ratio   = 0.71  (1)    : The ratio of NMDA to AMPA

    : Misc
    synapseID    = 0
    verboseLevel = 0
    LTPlasticity = 0
    rewiring     = 0
    synapseState = 2 : [0: Potential Synapse, 2: Functional Synapse]
}

COMMENT
The Verbatim block is needed to allow RNG.
ENDCOMMENT
VERBATIM
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM


ASSIGNED {
    : AMPAR assigned variables
    g_AMPA      (uS)
    i_AMPA      (nA)

    : NMDAR assigned variables
    g_NMDA      (uS)
    i_NMDA      (nA)
    mggate

    : Postsynaptic Ca2+ dynamics assigned variables
    : None

    : Release assigned variables
    Rstate      (1)  : recovered state {0=unrecovered, 1=recovered}
    tsyn        (ms)
    tsyn_fac    (ms) : the time of the last spike
    Psurv       (1)
    u           (1)  : running release probability (attention: u is event based based, so only valid at incoming events)

    : LTP/LTD assigned variables
    rng_GB
    val_rng_GB
    Theta_d_GB
    Theta_p_GB

    : Elimination assigned variables
    : None

    : Synaptogenesis assigned variables
    : None

    : Shared assigned variables
    v           (mV)
    vv          (mV)
    i           (nA)
    g           (uS)
    weight_AMPA
    weight_NMDA
    factor_AMPA
    factor_NMDA

    : Misc
    rng_rel

}

STATE {
    : AMPAR state variables to construct the dual-exponential profile
    A_AMPA : Decays with conductance tau_r_AMPA
    B_AMPA : Decays with conductance tau_d_AMPA

    : NMDAR state variables to construct the dual-exponential profile
    A_NMDA : Decays with conductance tau_r_NMDA
    B_NMDA : Decays with conductance tau_d_NMDA

    : Postsynaptic Ca2+ dynamics state variables
    cai_GB : Intracellular calcium concentration
    
    : LTP/LTD state variables
    rho_GB

    : Elimination state variables
    integrity_SE

    : Synaptogenesis state variables
    contact_SG
}


INITIAL {
    LOCAL tp_AMPA, tp_NMDA

    initialize_functional_synapse()

    : Time to peak of the conductances
    tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)
    tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA)

    : AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
    factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA)
    factor_AMPA = 1/factor_AMPA

    : NMDA Normalization factor - so that when t = tp_NMDA, gsyn = gpeak
    factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA)
    factor_NMDA = 1/factor_NMDA

    net_send(0, 1)
}


BEFORE BREAKPOINT {
    COMMENT
    In order to prevent integration errors, the function throws a single random
    number per time step and to ensure that the discontinuity in the rho SDE is
    handled properly.

    This implementation of the Graupner model is slightly different from the
    original. The choice of precomputing Theta_d_GB/Theta_p_GB might affect the
    convergence of the model.
    ENDCOMMENT

    val_rng_GB = nrand()
    Theta_d_GB = Theta(cai_GB - theta_d_GB)
    Theta_p_GB = Theta(cai_GB - theta_p_GB)
}


BREAKPOINT {
    SOLVE state METHOD cnexp

    : AMPAR 
    g_AMPA = gmax*(B_AMPA-A_AMPA)                                 : compute time varying conductance as the difference of state variables B_AMPA and A_AMPA
    i_AMPA = g_AMPA*(v-e)                                         : compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal

    :NMDAR
    mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (mg / 3.57 (mM))) : mggate kinetics - Jahr & Stevens 1990
    g_NMDA = gmax*(B_NMDA-A_NMDA) * mggate                        : compute time varying conductance as the difference of state variables B_NMDA and A_NMDA and mggate kinetics
    i_NMDA = g_NMDA*(v-e)                                         :compute the NMDA driving force based on the time varying conductance, membrane potential, and NMDA reversal

    g = g_AMPA + g_NMDA                                           : TODO Check why this line was here

    : Total current
    i = i_AMPA + i_NMDA

    vv = v
}


AFTER SOLVE {
    if( LTPlasticity == 1 ) {
        weight_AMPA = w0_GB + rho_GB*(w1_GB - w0_GB)
    } else {
        weight_AMPA = baseweight
    }
}


DERIVATIVE state{
    LOCAL rho0, a, b

    : AMPAR
    A_AMPA' = -A_AMPA/tau_r_AMPA
    B_AMPA' = -B_AMPA/tau_d_AMPA

    : NMDAR
    A_NMDA' = -A_NMDA/tau_r_NMDA
    B_NMDA' = -B_NMDA/tau_d_NMDA

    : Postsynaptic Ca2+ dynamics
    cai_GB'   = -cai_GB/tau_ca_GB

    : LTP/LTD
    :rho_GB' = ( -rho_GB*(1.0 - rho_GB)*(rho_star_GB - rho_GB)
    :            +gamma_p_GB*(1-rho_GB)*Theta_p_GB
    :            -gamma_d_GB*rho_GB*Theta_d_GB
    :            +Noise_GB() ) / (1000.0*tau_GB)
    rho0 = rho_GB
    a = ( -rho0*(1.0 - rho0)*(rho_star_GB - rho0)
          +gamma_p_GB*(1-rho0)*Theta_p_GB
          -gamma_d_GB*rho0*Theta_d_GB
          +Noise_GB() ) / (1000.0*tau_GB)
    b = ( -3.0*rho0*rho0 + 2*(1 + rho_star_GB)*rho0 -rho_star_GB
          -gamma_p_GB*Theta_p_GB
          -gamma_d_GB*Theta_d_GB ) / (1000.0*tau_GB)
    rho_GB' = a + b*(rho_GB - rho0)

    : Elimination
    :integrity_SE' = -integrity_SE/(tau_SE*(1+rho_GB))
    integrity_SE' = -integrity_SE/(tau_SE*(1+rho0))

    : Synaptogenesis
    contact_SG' = -contact_SG/tau_SG
}


NET_RECEIVE (weight){
    LOCAL result

    : Locals:
    : Psurv - survival probability of unrecovered state
    : tsyn - time since last surival evaluation.

    INITIAL{
        : TODO Check if this block is executed at the time of the first "true"
        :      spike or before the watch calls initialization. In the latter
        :      case, this is probably wrong.
    }

    if (flag == 1) {
        : Flag 1, Initialize watch calls
        WATCH (v > -30.0) 3
        WATCH (integrity_SE < theta_e_SE) 4
        WATCH (contact_SG > theta_g_SG) 5
        :printf("Flag 1, Initialize watch calls\n")

    } else if(flag == 2 && synapseState == 2){
        : Flag 2, presyn-induced calcium transient
        cai_GB = cai_GB + C_pre_GB
        :printf("Flag 2, presyn-induced calcium transient\n")

    } else if(flag == 3 && synapseState == 2){
        : Flag 3, postsyn-induced calcium transient
        cai_GB = cai_GB + C_post_GB
        :printf("Flag 3, postsyn-induced calcium transient\n")

    } else if(flag == 4 && rewiring == 1 && synapseState == 2){
        : Flag 4, eliminate synapse
        synapseState = 0
        initialize_potential_synapse()
        :printf("Flag 4, eliminate synapse\n")

    } else if(flag == 5 && rewiring == 1){
        : Flag 5, new functional synapse
        synapseState = 2
        initialize_functional_synapse()
        :printf("Flag 5, new functional synapse\n")

    } else if (flag == 0 && synapseState == 2){
        : Flag 0 (default) and synapse in the functional state
        weight_NMDA = baseweight * NMDA_ratio : Here for backward compatibility

        : calc u at event
        if (Fac > 0) {
            : Update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
            u = u*exp(-(t - tsyn_fac)/Fac)
            u = u + Use*(1-u)
        } else {
            u = Use
        }

        : tsyn_fac knows about all spikes, not only those that released
        : i.e. each spike can increase the u, regardless of recovered state.
        tsyn_fac = t

        : recovery
        if (Rstate == 0) {
            : Probability of survival of unrecovered state based on Poisson
            : recovery with rate 1/tau
            Psurv = exp(-(t-tsyn)/Dep)
            result = urand()
            if (result>Psurv) {
                : Recover
                Rstate = 1
                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf( "Recovered! %f at time %g: Psurv = %g, urand=%g\n",
                             synapseID, t, Psurv, result )
                    UNITSON
                }
            } else {
                : Survival must now be from this interval
                tsyn = t
                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf( "Failed to recover! %f at time %g: Psurv = %g, urand=%g\n", synapseID, t, Psurv, result )
                    UNITSON
                }
            }
        }

        if (Rstate == 1) {
            result = urand()
            if (result<u) {
                : release!
                tsyn = t
                Rstate = 0

                :printf("%.15f   %.15f\n", weight_AMPA, weight_NMDA)

                A_AMPA = A_AMPA + weight_AMPA*factor_AMPA
                B_AMPA = B_AMPA + weight_AMPA*factor_AMPA
                A_NMDA = A_NMDA + weight_NMDA*factor_NMDA
                B_NMDA = B_NMDA + weight_NMDA*factor_NMDA

                :printf("%.15f   %.15f\n", factor_AMPA, weight_AMPA*factor_AMPA)
                :printf("%.15f   %.15f\n", A_AMPA, B_AMPA)
                :printf("%.15f   %.15f\n", A_NMDA, B_NMDA)

                integrity_SE = integrity_SE + (1-integrity_SE)*factor_SE

                net_send(D_GB, 2)

                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf( "Release! %f at time %g: vals %g %g %g %g\n", synapseID, t, A_AMPA, weight_AMPA, factor_AMPA, weight )
                    UNITSON
                }
            } else {
                if( verboseLevel > 0 ) {
                    UNITSOFF
                    printf("Failure! %f at time %g: urand = %g\n", synapseID, t, result )
                    UNITSON
                }
            }
        }
    } else if (flag == 0 && synapseState == 0){
        : Flag 0 (default) and synapse in the potential state
        contact_SG = contact_SG + step_SG
    } else {
        : DO NOTHING
    }
}


PROCEDURE initialize_functional_synapse() {
    Rstate=1
    tsyn_fac=0
    u=u0

    A_AMPA = 0
    B_AMPA = 0

    A_NMDA = 0
    B_NMDA = 0

    : Postsynaptic Ca2+ concentration
    cai_GB = 0

    : LTP/LTD
    rho_GB = (baseweight - w0_GB) / (w1_GB - w0_GB)

    : Elimination
    integrity_SE = 1

    contact_SG = 0
}


PROCEDURE initialize_potential_synapse() {
    Rstate=1
    tsyn_fac=0
    u=u0

    A_AMPA = 0
    B_AMPA = 0

    A_NMDA = 0
    B_NMDA = 0

    : Postsynaptic Ca2+ concentration
    cai_GB = 0

    : LTP/LTD
    rho_GB = 0

    : Elimination
    integrity_SE = 0

    contact_SG = 0
}


FUNCTION Theta(x) {
    if (x < 0.0) {
        Theta = 0.0
    } else {
        Theta = 1.0
    }
}


FUNCTION Noise_GB() {
    Noise_GB = sigma_GB*sqrt(tau_GB)*sqrt(Theta_p_GB + Theta_d_GB)*val_rng_GB

    :printf("t = %f    val_rng = %f    Noise = %f    Theta_d = %f    Theta_p = %f\n", t, val_rng_GB, Noise_GB, Theta_d_GB, Theta_p_GB)
}


PROCEDURE setRNG() {
    VERBATIM
    /**
     * This function takes a NEURON Random object declared in hoc and makes it
     * usable by this mod file.
     */
    void** pv1 = (void**)(&_p_rng_rel);
    void** pv2 = (void**)(&_p_rng_GB);
    if( ifarg(2)) {
        *pv1 = nrn_random_arg(1);
        *pv2 = nrn_random_arg(2);
    } else {
        *pv1 = (void*)0;
        *pv2 = (void*)0;
    }
    ENDVERBATIM
}


FUNCTION urand() {
    VERBATIM
    double value;
    if (_p_rng_rel) {
        /*
        :Supports separate independent but reproducible streams for
        : each instance. However, the corresponding hoc Random
        : distribution MUST be set to Random.uniform(0,1)
        */
        value = nrn_random_pick(_p_rng_rel);
        //printf("random stream for this simulation = %lf\n",value);
        return value;
    }else{
    ENDVERBATIM
        : the old standby. Cannot use if reproducible parallel sim
        : independent of nhost or which host this instance is on
        : is desired, since each instance on this cpu draws from
        : the same stream
        value = scop_random(1)
    VERBATIM
    }
    ENDVERBATIM
    urand = value
}


FUNCTION nrand() {
    VERBATIM
    double value;
    if (_p_rng_GB) {
        /*
        :Supports separate independent but reproducible streams for
        : each instance. However, the corresponding hoc Random
        : distribution MUST be set to Random.normal(0,1)
        */
        value = nrn_random_pick(_p_rng_GB);
        //printf("random stream for this simulation = %lf\n",value);
        return value;
    }else{
    ENDVERBATIM
        : the old standby. Cannot use if reproducible parallel sim
        : independent of nhost or which host this instance is on
        : is desired, since each instance on this cpu draws from
        : the same stream
        value = scop_random(1)
    VERBATIM
    }
    ENDVERBATIM
    nrand = value
}


FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}


FUNCTION toggleLTPlasticity() {
    LTPlasticity = 1-LTPlasticity
}


FUNCTION toggleRewiring() {
    rewiring = 1-rewiring
}


FUNCTION bbsavestate() {
        bbsavestate = 0
VERBATIM
#ifdef ENABLE_SAVE_STATE
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng_rel) {
                // tell how many items need saving
                if (*xdir == -1. ) { *xdir = 2.0; return 0.0; }

                // save the value(s)
                else if (*xdir == 0.) {
                        xval[0] = (double) nrn_get_random_sequence(_p_rng_GB);
                        xval[1] = (double) nrn_get_random_sequence(_p_rng_rel);
                } else{  //restore the value(s)
                        nrn_set_random_sequence(_p_rng_GB, (long)(xval[0]));
                        nrn_set_random_sequence(_p_rng_rel, (long)(xval[1]));
                }
        }

        //if( synapseID == 104211 ) { verboseLevel = 1; }
#endif
ENDVERBATIM
}

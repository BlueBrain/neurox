COMMENT
/**
 * @file StdpWADoublet.mod
 * @brief 
 * @author king
 * @date 2011-06-22
 * @remark Copyright © BBP/EPFL 2005-2011; All rights reserved. Do not distribute without further notice.
 */
ENDCOMMENT

COMMENT
Spike Timing Dependent Weight Adjuster
Matt Perich, EPFL, 2011

NOTE: This file has some hooks built in for implementing
 the Ca2+ rule Eilif and I worked on. If the Ca rule isn't
 configured, then the params have no effect on how the rule behaves.
ENDCOMMENT

NEURON {
    POINT_PROCESS StdpWADoublet
    RANGE wmax, wmin, aPlus, aMinus, tauPlus, tauMinus, on
    RANGE allow_update_on_post, verbose
    RANGE r, o
    RANGE m, sigmed, aPlusNew
    RANGE weight, synapseID
    POINTER wsyn, eta
}

ASSIGNED {
    wsyn		    	: weight of the synapse
    eta
    aPlusNew
    weight
}

STATE {
    r   : tracks postsynaptic spikes with a decaying exponential
    o   : tracks presynaptic spikes with a decaying exponential
}

INITIAL {
    r = 0
    o = 0
    net_send(0,1)
}

PARAMETER {
    tauPlus  = 20	(ms) : decay time for LTP part ( values from           )
    tauMinus = 20	(ms) : decay time for LTD part ( Song and Abbott, 2001 )
    wmax     = 5.0            : min and max values
    wmin     = 0              :   of synaptic weight
    aPlus    = 0.01           : amplitude of LTP steps
    aMinus   = 0.01	     : amplitude of LTD steps
    m	     = 0.08	     : slope of sigmoidal function for aLTP changing
    sigmed   = 0.5 	     : inflection point of sigmoid
    on       = 1		     : allows learning to be turned on and off
    verbose  = 0              : if true, will write info to the screen
    allow_update_on_post = 1 : if this is true, we update on receiving both pre-/post-synaptic spikes
                             :   if false weight updates accumulate and get set only for pre-synaptic
    synapseID = 0
}

BREAKPOINT {
    SOLVE state METHOD cnexp
}

DERIVATIVE state {
    r' = -r/tauPlus
    o' = -o/tauMinus
}

NET_RECEIVE (w) {

    if(flag==1) {           : initialize the weight variable
    
        weight = wsyn
	aPlusNew = aPlus
    
    } else {

        :aPlusNew = aPlus/(1.0+exp((sigmed-eta)/m))
    
    	if (w >= 0) {       : this is a pre-synaptic spike

   	    o = o + 1

            if (verbose == 1) {
                printf("D: weight=%0.4f, change=%0.5f, r=%0.4f, o=%0.3f\n",weight,aMinus*(weight-wmin)*r,r,o)
            }

            weight = weight - aMinus * (weight - wmin) * r
            
    	} else {		    : this is a post-synaptic spike

    	    r = r + 1

            if (verbose == 1) {
                printf("P: weight=%0.4f, change=%0.5f, o=%0.4f, r=%0.3f\n",weight,aPlusNew*(wmax-weight)*o,o,r)
            }
	
            weight = weight + aPlusNew * (wmax - weight) * o

    	}
    	
    	if (on) { : turning plasticity on or off, then clipping to max/min
    	    if (weight > wmax) {
    	        weight = wmax
    	    }
            if (weight < wmin) {
                weight = wmin
    	    }
    		    
    	    if (w >= 0 || allow_update_on_post) {
    	        wsyn = weight : adjust the weight pointer to the synapse
    	    }
    		
    	}
    }
}

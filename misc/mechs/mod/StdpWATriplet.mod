COMMENT
/**
 * @file StdpWATriplet.mod
 * @brief 
 * @author king
 * @date 2011-06-22
 * @remark Copyright © BBP/EPFL 2005-2011; All rights reserved. Do not distribute without further notice.
 */
ENDCOMMENT


COMMENT
Spike Timing Dependent Weight Adjuster
Uses a triplet rule as proposed by Pfister/Gerstner
Matt Perich, EPFL, 2011
ENDCOMMENT

NEURON {
    POINT_PROCESS StdpWATriplet
    RANGE wmin, wmax, a2Plus, a3Plus, a2Minus, a3Minus
    RANGE tauPlus, tauMinus, tauX, tauY
    RANGE allow_update_on_post, on, synapseID
    RANGE r1, r2, o1, o2
    RANGE weight
    POINTER wsyn
    POINTER eta : For calcium rule... not used currently
}

ASSIGNED {
    wsyn : weight of the synapse
    weight
    eta
}

STATE {
    r1   : descriptions
    r2   : 
    o1   : 
    o2   : 
}

INITIAL {
    r1 = 0
    r2 = 0
    o1 = 0
    o2 = 0
    net_send(0,1)
}

PARAMETER {
    wmin = 0.0
    wmax = 1.0
    tauPlus  = 16.8	(ms)     : decay time for LTP part ( values from           )
    tauMinus  = 33.7	(ms)     : decay time for LTD part ( Song and Abbott, 2001 )
    tauX    = 101
    tauY    = 125
    a2Plus  = 0.0000000000005	    	 : amplitude of LTP steps
    a3Plus  = 0.0000062
    a2Minus = 0.000007	     : amplitude of LTD steps
    a3Minus = 0.00000023
    synapseID = 0
    on      = 1		         : allows learning to be turned on and off
    allow_update_on_post = 1 : if this is true, we update on receiving both pre-/post-synaptic spikes
                             :   if false weight updates are accumulated and applied only for pre-synaptic
}

BREAKPOINT {
    SOLVE state METHOD cnexp
}

DERIVATIVE state {
    r1' = -r1/tauPlus
    r2' = -r2/tauX
    o1' = -o1/tauMinus
    o2' = -o2/tauY
}

NET_RECEIVE (w) {

    if(flag==1) {           : initialize the weight variable
    
        weight = wsyn
    
    } else {

    	if (w >= 0) {       : this is a pre-synaptic spike
            weight = weight - o1 * (a2Minus + a3Minus * r2)
    	    r1 = r1 + 1.0
            r2 = r2 + 1.0
    	} else {		    : this is a post-synaptic spike	
            weight = weight + r1 * (a2Plus + a3Plus * o2)
    	    o1 = o1 + 1.0
            o2 = o2 + 1.0
    	}
    	
    	if (on) { 

            if (weight < wmin) {
                weight = wmin
    	    }
    		    
    	    if (w >= 0 || allow_update_on_post) {
    	        wsyn = weight
    	    }
    		
    	}
    }
}

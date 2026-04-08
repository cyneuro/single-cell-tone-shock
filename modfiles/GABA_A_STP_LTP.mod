COMMENT
/**
 * @file GABA_A_STP_LTP.mod
 * @brief GABAA receptor with presynaptic short-term plasticity and long-term plasticity
 * @author Modified from GABA_A_STP.mod and AMPA_NMDA_STP_LTP.mod
 * @date 2025-10-02
 * @remark Combines GABAA receptor conductance with STP and calcium-dependent LTP/LTD
 */
ENDCOMMENT

TITLE GABAA receptor with presynaptic short-term plasticity and long-term plasticity

COMMENT
GABAA receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al. 2002, deterministic version.
Long-term plasticity based on Shouval et al. 2002a, 2002b with calcium dynamics
ENDCOMMENT

NEURON {
    THREADSAFE

    POINT_PROCESS GABA_A_STP_LTP
    USEION ca READ eca, ica
    RANGE initW     : synaptic scaler for large scale networks
    RANGE tau_r_GABAA, tau_d_GABAA
    RANGE Use, Dep, Fac, u0
    RANGE gmax, gmax_GABAA
    RANGE i, g, e_GABAA
    NONSPECIFIC_CURRENT i
    RANGE synapseID, verboseLevel
    RANGE conductance
    GLOBAL nc_type_param
    RANGE record_use, record_Pr
    
    : Long-term plasticity parameters
    RANGE Cainf, tauCa, Icatotal, volume_CR
    RANGE ICag, P0g, fCag, gamma_ca_CR, tau_effca_GB
    RANGE lambda1, lambda2, threshold1, threshold2
    RANGE fmax, fmin, Wmax, Wmin, limitW
    RANGE W, cai_CR, effcai_GB, LTP_on
    RANGE dep_GB, pot_GB
}

PARAMETER {
    initW        = 1.0        : synaptic weight scaler
    tau_r_GABAA  = 0.2   (ms) : dual-exponential conductance profile
    tau_d_GABAA  = 8     (ms) : IMPORTANT: tau_r < tau_d
    Use          = 1.0   (1)  : Utilization of synaptic efficacy
    Dep          = 100   (ms) : relaxation time constant from depression
    Fac          = 10    (ms) : relaxation time constant from facilitation
    e_GABAA      = -75   (mV) : GABAA reversal potential
    gmax         = .001  (uS) : weight conversion factor (from nS to uS)
    u0           = 0          : initial value of u, release probability
    synapseID    = 0
    verboseLevel = 0
    conductance  = 0.0
    nc_type_param = 7
    
    : Calcium dynamics for LTP (matching AMPA system)
    Cainf = 70e-6 (mM)        : Baseline calcium concentration
    tauCa = 12 (ms)           : Fast calcium removal time constant  
    gamma_ca_CR = 0.04 (1)    : Percent of free calcium (not buffered)
    tau_effca_GB = 200 (ms)   : Slow effective calcium time constant
    volume_CR = 0.087 (um3)   : Spine volume (matching AMPA)
    
    P0g = .01
    fCag = .024
    
    : LTP parameters
    lambda1 = 15.0            : LTP strength parameter (potentiation rate)
    lambda2 = 0.001            : LTD strength parameter (depression rate)
    threshold1 = 2 (us/liter) : LTD threshold (adjusted to match AMPA scale)
    threshold2 = 3 (us/liter) : LTP threshold (adjusted to match AMPA scale)
    LTP_on = 1                 : LTP switch (1=on, 0=off)
    
    fmax = 3.0            : maximum weight factor
    fmin = 0.8            : minimum weight factor
    
    k = 0.01              : calcium current scaling factor
}

UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
    (uS) = (microsiemens)
    (mM) = (milli/liter)
    (uM) = (micro/liter)
    (us/liter) = (micro/liter)
    FARADAY = 96485 (coul)
    pilocal = 3.141592 (1)
}

ASSIGNED {
    v (mV)
    eca (mV)
    ica (nA)
    
    i (nA)
    g (uS)
    gmax_GABAA (uS)
    factor_GABAA
    record_use
    record_Pr
    
    limitW
    Wmax
    Wmin
    
    ICag (nA)
    Icatotal (nA)
    dep_GB (1)
    pot_GB (1)
}

STATE { 
    A_GABAA       : GABAA state variable - decays with tau_r_GABAA
    B_GABAA       : GABAA state variable - decays with tau_d_GABAA
    cai_CR        : intracellular calcium concentration (mM) - fast calcium dynamics (tau=12ms)
    effcai_GB     : effective calcium for plasticity (us/liter) - slow integration of calcium elevations above baseline (tau=200ms) used for LTP/LTD threshold detection <1e-3>
    W             : synaptic weight
}

INITIAL{
    LOCAL tp_GABAA

    A_GABAA = 0
    B_GABAA = 0
    cai_CR = Cainf
    effcai_GB = 0
    W = initW
    if (LTP_on == 0) {
        W = initW
    }

    tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA)
    factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA)
    factor_GABAA = 1/factor_GABAA

    gmax_GABAA = initW * gmax

    Wmax = fmax*initW
    Wmin = fmin*initW
    limitW = 1

    record_use = u0
    record_Pr = u0
    dep_GB = 0
    pot_GB = 0
    net_send(0, 1)
}

BREAKPOINT {
    : Limit weight changes based on effective calcium (matching AMPA system)
    if ((eta((effcai_GB / 1000))*(lambda1*omega((effcai_GB / 1000), threshold1, threshold2)-lambda2*W))>0 && W>=Wmax) {
        limitW = 1e-12
    } else if ((eta((effcai_GB / 1000))*(lambda1*omega((effcai_GB / 1000), threshold1, threshold2)-lambda2*W))<0 && W<=Wmin) {
        limitW = 1e-12
    } else {
        limitW = 1
    }
    
    SOLVE release METHOD cnexp
    
    g = gmax_GABAA*(B_GABAA-A_GABAA)
    i = W*g*(v - e_GABAA)
    
    ICag = P0g*g*(v - eca)
    Icatotal = ICag + k*ica*4*pilocal*((15/2)^2)*(0.01)
    
    : Clamp W to limits
    if (W > Wmax) {
        W = Wmax
    } else if (W < Wmin) {
        W = Wmin
    }
}

DERIVATIVE release {
    : Weight dynamics with LTP/LTD (using effective calcium like AMPA)
    if (LTP_on == 1) {
        W' = limitW*eta((effcai_GB / 1000))*(lambda1*omega((effcai_GB / 1000), threshold1, threshold2)-lambda2*W)
    } else {
        W' = 0
    }
    
    : GABAA receptor kinetics
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
    
    : Two-stage calcium dynamics (matching AMPA system)
    : Fast intracellular calcium: increases from GABA-mediated calcium influx, decays to baseline (12ms)
    cai_CR' = - (1e-9)*Icatotal*gamma_ca_CR/((1e-15)*volume_CR*2*FARADAY) - (cai_CR - Cainf)/tauCa
    : Slow effective calcium: integrates calcium elevations above baseline for plasticity decisions (200ms)
    effcai_GB' = - effcai_GB/tau_effca_GB + (cai_CR - Cainf)
}

NET_RECEIVE (weight, R, u, tsyn (ms)){
    LOCAL Pr, weight_GABAA
    INITIAL{
        R = 1
        u = u0
        tsyn = t
    }

    if (flag == 0) {
        : Short-term plasticity calculations
        if (Fac > 0) {
            u = u*exp(-(t - tsyn)/Fac)
            u = u + Use*(1-u)
        } else {
            u = Use
        }

        if (Dep > 0) {
            R = 1 - (1-R) * exp(-(t-tsyn)/Dep)
            Pr = u * R
            R = R - u * R
        } else {
            Pr = u 
        }

        record_use = u
        record_Pr = Pr

        if( verboseLevel > 0 ) {
            printf("Synapse %f at time %g: R = %g Pr = %g\n", synapseID, t, R, Pr )
        }

        tsyn = t

        weight_GABAA = Pr*weight*factor_GABAA
        A_GABAA = A_GABAA + weight_GABAA
        B_GABAA = B_GABAA + weight_GABAA

        if( verboseLevel > 0 ) {
            printf( " vals %g %g %g %g\n", A_GABAA, weight_GABAA, factor_GABAA, weight )
        }
    } else if (flag == 1) {
        WATCH (effcai_GB > threshold1) 2
        WATCH (effcai_GB < threshold1) 3
        WATCH (effcai_GB > threshold2) 4
        WATCH (effcai_GB < threshold2) 5
    } else if (flag == 2) {
        dep_GB = 1
    } else if (flag == 3) {
        dep_GB = 0
    } else if (flag == 4) {
        pot_GB = 1
    } else if (flag == 5) {
        pot_GB = 0
    }
}

: LTP/LTD functions adapted from Shouval et al. 2002
FUNCTION eta(Cani (us/liter)) {
    LOCAL taulearn, P1, P2, P4, Cacon
    P1 = 0.1
    P2 = P1*1e-4
    P4 = 1
    Cacon = Cani*1e3
    taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4
    eta = 1/taulearn*0.001
}

FUNCTION omega(Cani (us/liter), threshold1 (us/liter), threshold2 (us/liter)) {
    LOCAL r, mid, Cacon
    Cacon = Cani*1e3
    r = (threshold2-threshold1)/2
    mid = (threshold1+threshold2)/2
    if (Cacon <= threshold1) { 
        omega = 0
    } else if (Cacon >= threshold2) {
        omega = 1/(1+50*exp(-50*(Cacon-threshold2)))
    } else {
        omega = -sqrt(r*r-(Cacon-mid)*(Cacon-mid))
    }
}

FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}
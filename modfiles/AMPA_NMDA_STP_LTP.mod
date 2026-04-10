COMMENT
/**
 * @file LTP_simple.mod
 * @brief AMPA and NMDA receptor with LTP based on calcium-dependent plasticity rules
 * @remark Adapted from AMPA_NMDA_STP.mod and pyrD2pyrD_STFD_new.mod
 */
ENDCOMMENT


TITLE AMPA and NMDA receptor with LTP


COMMENT
AMPA and NMDA receptor conductance using a dual-exponential profile
with calcium-dependent long-term potentiation.
ENDCOMMENT


NEURON {
    THREADSAFE

    POINT_PROCESS AMPA_NMDA_STP_LTP
    : AMPA/NMDA dual-exponential conductance time constants (rise and decay)
    RANGE tau_r_AMPA, tau_d_AMPA, tau_r_NMDA, tau_d_NMDA
    : Short-term plasticity parameters: utilization, depression timescale, facilitation timescale, initial utilization
    RANGE Use, Dep, Fac, u0
    : Synaptic currents (total, AMPA, NMDA) and conductances + reversal potential
    RANGE i, i_AMPA, i_NMDA, g_AMPA, g_NMDA, g, e
    : Maximum conductances for AMPA and NMDA receptors
    RANGE gmax_AMPA, gmax_NMDA
    NONSPECIFIC_CURRENT i
    : Utility parameters: synapse identifier and debug verbosity level
    RANGE synapseID, verboseLevel
    : For LTP
    USEION ca READ eca
    : LTP weight parameters: initial weight, current weight, weight bounds, learning rate parameters, calcium thresholds, LTP on/off flag
    RANGE initW, W, Wmax, Wmin, lambda1, lambda2, threshold1, threshold2, LTP_on
    : Calcium dynamics parameters: resting calcium, pool diameter, valence, calcium time constant, NMDA current fraction, calcium buffering, calcium accumulation factor, NMDA calcium current, total calcium current, local calcium concentration
    RANGE Cainf, pooldiam, z, tauCa, P0, fCa, Afactor, ICa, iCatotal, capoolcon
    : Recorded variables for debugging: running utilization and release probability
    RANGE record_use, record_Pr
}


PARAMETER {
    : AMPA kinetics: rise time for dual-exponential conductance profile
    tau_r_AMPA = 0.2   (ms)
    : AMPA kinetics: decay time for dual-exponential conductance profile (IMPORTANT: tau_r < tau_d)
    tau_d_AMPA = 1.7   (ms)
    : NMDA kinetics: rise time for dual-exponential conductance profile
    tau_r_NMDA = 0.29  (ms)
    : NMDA kinetics: decay time for dual-exponential conductance profile (IMPORTANT: tau_r < tau_d)
    tau_d_NMDA = 43    (ms)
    : Short-term plasticity: initial utilization of synaptic efficacy (fraction of available resources released per spike)
    Use = 1.0          (1)
    : Short-term plasticity: time constant for recovery from synaptic depression
    Dep = 100          (ms)
    : Short-term plasticity: time constant for decay of facilitation
    Fac = 10           (ms)
    : Reversal potential for both AMPA and NMDA receptors
    e = 0              (mV)
    : Maximum NMDA conductance; scaling factor converts from nanosiemens to microsiemens
    gmax_NMDA = .001   (uS)
    : Maximum AMPA conductance; scaling factor converts from nanosiemens to microsiemens
    gmax_AMPA = .001   (uS)
    : Relative strength of NMDA current compared to AMPA (NMDA conductance = AMPA conductance * NMDA_ratio)
    NMDA_ratio = 0.71  (1)
    : Initial value of u (utilization running variable) at rest, determines baseline release probability
    u0 = 0             (1)
    : Unique identifier for this synapse instance (for logging/debugging)
    synapseID = 0      (1)
    : Debug flag: enables detailed printf output when > 0
    verboseLevel = 0   (1)

    : ===== LTP PARAMETERS =====
    : Initial synaptic weight (before any plasticity)
    initW = 1.0        (1)
    : LTP learning rate: magnitude of weight increase when calcium is elevated (scales potentiation)
    lambda1 = 40       (1)
    : LTP learning rate: magnitude of weight decrease when calcium is low (scales depression)
    lambda2 = 0.03     (1)
    : Lower calcium threshold: below this level depotentiation occurs
    threshold1 = 0.4   (uM)
    : Upper calcium threshold: above this level potentiation occurs
    threshold2 = 0.55  (uM)
    : Maximum weight bound: weight multiplied by this factor to set upper limit (sets Wmax = fmax * initW)
    fmax = 3           (1)
    : Minimum weight bound: weight multiplied by this factor to set lower limit (sets Wmin = fmin * initW)
    fmin = 0.8         (1)
    : Resting (baseline) intracellular calcium concentration
    Cainf = 50e-6      (mM)
    : Diameter of calcium pool compartment (used to calculate volume for calcium concentration changes)
    pooldiam = 1.8172  (micrometer)
    : Valence of calcium ion (charge)
    z = 2              (1)
    : Time constant for calcium clearance/decay back to resting level
    tauCa = 50         (ms)
    : Fraction of NMDA current that is carried by calcium (P0 * g_NMDA * driving force = ICa)
    P0 = 0.015         (1)
    : Calcium accumulation factor: scales NMDA calcium current contribution to calcium pool
    fCa = 0.024        (1)
    : ===== LTP CONTROL FLAG =====
    : Enable/disable LTP weight changes: 1 = LTP active, 0 = LTP inactive (weight frozen)
    LTP_on = 1         (1)
}


UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
    (uS) = (microsiemens)
    FARADAY = 96485 (coul)
    pilocal = 3.141592 (1)
}


ASSIGNED {
    v (mV)
    i (nA)
    i_AMPA (nA)
    i_NMDA (nA)
    g_AMPA (uS)
    g_NMDA (uS)
    g (uS)
    factor_AMPA
    factor_NMDA
    eca (mV)
    ICa (mA)
    Afactor (mM/ms/nA)
    iCatotal (mA)
    Wmax
    Wmin
    limitW
    record_use
    record_Pr
}


STATE {
    A_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_r_AMPA
    B_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_d_AMPA
    A_NMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_r_NMDA
    B_NMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_d_NMDA
    capoolcon    : local calcium concentration
    W            : synaptic weight
}


INITIAL{
    LOCAL tp_AMPA, tp_NMDA

    A_AMPA = 0
    B_AMPA = 0

    A_NMDA = 0
    B_NMDA = 0

    tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA) :time to peak of the conductance
    tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA) :time to peak of the conductance

    factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) :AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
    factor_AMPA = 1/factor_AMPA

    factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA) :NMDA Normalization factor - so that when t = tp_NMDA, gsyn = gpeak
    factor_NMDA = 1/factor_NMDA

    capoolcon = Cainf
    W = initW
    Wmax = fmax * initW
    Wmin = fmin * initW
    Afactor = 1/(z*FARADAY*4/3*pilocal*(pooldiam/2)^3)*(1e6)
    limitW = 1

    record_use = u0
    record_Pr = u0
}


BREAKPOINT {
    if ((eta(capoolcon)*(lambda1*omega(capoolcon, threshold1, threshold2)-lambda2*W))>0&&W>=Wmax) {
        limitW=1e-12
    } else if ((eta(capoolcon)*(lambda1*omega(capoolcon, threshold1, threshold2)-lambda2*W))<0&&W<=Wmin) {
        limitW=1e-12
    } else {
        limitW=1
    }

    SOLVE state METHOD cnexp
    g_AMPA = gmax_AMPA*(B_AMPA-A_AMPA) :compute time varying conductance as the difference of state variables B_AMPA and A_AMPA
    g_NMDA = gmax_NMDA*(B_NMDA-A_NMDA) * sfunc(v) :compute time varying conductance as the difference of state variables B_NMDA and A_NMDA and mggate kinetics
    g = g_AMPA + g_NMDA
    i_AMPA = g_AMPA*(v-e) :compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal
    i_NMDA = g_NMDA*(v-e) :compute the NMDA driving force based on the time varying conductance, membrane potential, and NMDA reversal
    i = (i_AMPA + i_NMDA) * W

    ICa = P0 * g_NMDA * (v - eca)
}


DERIVATIVE state{
    A_AMPA' = -A_AMPA/tau_r_AMPA
    B_AMPA' = -B_AMPA/tau_d_AMPA
    A_NMDA' = -A_NMDA/tau_r_NMDA
    B_NMDA' = -B_NMDA/tau_d_NMDA
    capoolcon' = -fCa * Afactor * ICa + (Cainf - capoolcon) / tauCa
    W' = limitW * eta(capoolcon) * (lambda1 * omega(capoolcon, threshold1, threshold2) - lambda2 * W) * LTP_on
}


NET_RECEIVE (weight, weight_AMPA, weight_NMDA, R, Pr, u, tsyn (ms)){
    weight_AMPA = weight
    weight_NMDA = weight * NMDA_ratio

    INITIAL{
        R=1
        u=u0
        tsyn=t
    }

    : flag == 0, i.e. a spike has arrived

    : calc u at event-
    if (Fac > 0) {
        u = u*exp(-(t - tsyn)/Fac) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
        u = u + Use*(1-u) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    } else {
        u = Use
    }

    if (Dep > 0) {
        R  = 1 - (1-R) * exp(-(t-tsyn)/Dep) :Probability R for a vesicle to be available for release, analogous to the pool of synaptic
                                        :resources available for release in the deterministic model. Eq. 3 in Fuhrmann et al.
        Pr = u * R                      :Pr is calculated as R * u (running value of Use)
        R  = R - u * R                  :update R as per Eq. 3 in Fuhrmann et al.
    } else {
        Pr = u 
    }

    record_use = u
    record_Pr = Pr

    if( verboseLevel > 0 ) {
        printf("Synapse %f at time %g: R = %g Pr = %g\n", synapseID, t, R, Pr )
    }

    tsyn = t

    A_AMPA = A_AMPA + Pr*weight_AMPA*factor_AMPA
    B_AMPA = B_AMPA + Pr*weight_AMPA*factor_AMPA
    A_NMDA = A_NMDA + Pr*weight_NMDA*factor_NMDA
    B_NMDA = B_NMDA + Pr*weight_NMDA*factor_NMDA

    if( verboseLevel > 0 ) {
        printf( " vals %g %g %g %g\n", A_AMPA, weight_AMPA, factor_AMPA, weight )
    }
}


FUNCTION sfunc (v (mV)) {
    UNITSOFF
    sfunc = 1/(1+0.33*exp(-0.06*v))
    UNITSON
}


FUNCTION eta(Cani (mM)) {
    LOCAL taulearn, P1, P2, P4, Cacon
    P1 = 0.1
    P2 = P1*1e-4
    P4 = 1
    Cacon = Cani*1e3
    taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4
    eta = 1/taulearn*0.001
}


FUNCTION omega(Cani (mM), threshold1 (uM), threshold2 (uM)) {
    LOCAL r, mid, Cacon
    Cacon = Cani*1e3
    r = (threshold2-threshold1)/2
    mid = (threshold1+threshold2)/2
    if (Cacon <= threshold1) { omega = 0}
    else if (Cacon >= threshold2) { omega = 1/(1+50*exp(-50*(Cacon-threshold2)))}
    else {omega = -sqrt(r*r-(Cacon-mid)*(Cacon-mid))}
}


FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}
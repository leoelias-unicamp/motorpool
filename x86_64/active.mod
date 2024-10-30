TITLE active.mod
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX active
        USEION ca READ eca WRITE ica
	NONSPECIFIC_CURRENT il
        RANGE gcabar, gl, el, gca, theta
        GLOBAL pinf, ptau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gcabar = .000038 (S/cm2)	<0,1e9>
	
	gl = .0003 (S/cm2)	<0,1e9>
        el = 0.0 (mV)

	theta = 25.0 (mV)

}
 
STATE {
        p
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        eca (mV)

	gca (S/cm2)
        ica (mA/cm2)
        il (mA/cm2)
        pinf
	ptau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = gcabar*p
	ica = gca*(v - eca)
        il = gl*(v - el)
}
  
INITIAL {
	rates(v)
	p = pinf
}

? states
DERIVATIVE states {  
        rates(v)
        p' =  (pinf-p)/ptau
}

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL theta_p, k_p
        TABLE pinf, ptau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
                :"p" sodium activation system
	theta_p = theta
	k_p = -7.0
	ptau = 20.0
        pinf = 1/(1 + exp((v-theta_p)/k_p))
}

UNITSON

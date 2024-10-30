TITLE traub.mod   modified traub model
 
COMMENT
 This is a modification of the Traub neuron model, which includes three ionic
channels: sodium, fast potassium, and slow potassium. The parameters were set
so as to reproduce the behavior of a spinal motoneuron. 
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX traub
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
        RANGE gnabar, gkfbar, gksbar, gl, el, gna, gkf, gks, vtraub
        GLOBAL minf, hinf, ninf, rinf, mtau, htau, ntau, rtau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .030 (S/cm2)	<0,1e9>
        gkfbar = .004 (S/cm2)	<0,1e9>
	gksbar = .016 (S/cm2)	<0,1e9>
 	
	gl = .0003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)

	vtraub = 0.0 (mV)
}
 
STATE {
        m h n r
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gkf (S/cm2)
	gks (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf rinf
	mtau (ms) htau (ms) ntau (ms) rtau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
        gkf = gkfbar*n*n*n*n
	gks = gksbar*r*r
	ik = (gkf+gks)*(v-ek)
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	r = rinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
	r' = (rinf-r)/rtau
}

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum, v2
        TABLE minf, mtau, hinf, htau, ninf, ntau, rinf, rtau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF

	v2 = v - vtraub

                :"m" sodium activation system
        alpha = .32 * vtrap(13-v2,4)
        beta = .28 * vtrap(v2-40,5)
        sum = alpha + beta
	mtau = 1/sum
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .128 * exp((17-v2)/18)
        beta = 1 / (exp((40-v2)/5) + 1)
        sum = alpha + beta
	htau = 1/sum
        hinf = alpha/sum
                :"n" fast potassium activation system
        alpha = .032*vtrap(15-v2,5) 
        beta = .5*exp((10-v2)/40)
	sum = alpha + beta
        ntau = 1/sum
        ninf = alpha/sum
		:"r" slow potassium activation system
	alpha = 3.5/exp(((55-v2)/4) + 1)
	beta = 0.025
	sum = alpha + beta
	rtau = 1/sum
	rinf = alpha/sum 
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

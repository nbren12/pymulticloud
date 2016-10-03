"""

Fortran Code
============================================================

READ(69,*) qbar
READ(69,*) dbar
READ(69,*) pbar
READ(69,*) qr02
READ(69,*) tau_e
READ(69,*) m0
READ(69,*) cape0
READ(69,*) capebar
READ(69,*) fceq
READ(69,*) fdeq
READ(69,*) fseq
READ(69,*) rstoch


READ(69,*) tau01
READ(69,*) tau02
READ(69,*) tau10
READ(69,*) tau12
READ(69,*) tau20
READ(69,*) tau30
READ(69,*) tau23

READ(69,*) alpha_bar
READ(69,*) l
READ(69,*) c
READ(69,*) moist0
"""
import numpy as np
from sympy import *
from scipy.optimize import fsolve

hour_s = 3600.0
day_s  = 24 * hour_s
alpha_bar = 15.3061
Te = 8.33 * hour_s
c  = 50.0
Le = c * Te
Hm  = 5000.0 / Le # this is weird, why non dim by horiz length scale?
Ht  = 16e3   / Le
Hb  = 500 / Le


#######################################################################
#                     Use sympy to write the RHS                      #
#######################################################################
tstpi = 2 * sqrt(2) /pi
alps, alpc, qbar, tc, xis, xic, alpha2, temteb_bar=symbols('alpha_s alpha_c qbar tau_conv xis xic alpha3 theta_eb_m_theta_em')# in code alpha3 is used in theta_em, in paper it is called alpha3
qr01, sst_bar, T0=symbols('qr01 theta_ebs_m_theta_eb moist0')

# Stochastic Time Scales
tau01, tau02, tau10, tau12, tau20, tau30, tau23 =   symbols('tau01 tau02 tau10 tau12 tau20 tau30 tau23')

m0, mu, tau_e, capebar, qr02=symbols('m0 mu tau_e capebar qr02')

sigS, sigD, sigC = symbols('sigma_s sigma_d sigma_c')
fseq, fdeq, fceq = symbols('sigma_s_bar sigma_d_bar sigma_c_bar')
u1, u2, Hs, teb, tebs, t1, t2, q = symbols('u1 u2 Hs theta_eb thteb_st theta1 theta2 q')


Rad1 = sympify('qr01 + theta1 / tau_r')
Rad2 = sympify('qr02 + theta2 / tau_r')

CAPE = sympify('capebar + rstoch *( theta_eb - 1.7 * (theta1 + alpha2*theta2))')
CAPEl = sympify('capebar + rstoch *( theta_eb - 1.7 * (theta1 + alpha4*theta2))')

CAPEm = sympify('a1 * theta_eb + a2 *q -a0 *( theta1 +alpha2 * theta2 )')
Hd = sigD*qbar + (sigD / fdeq / tc) * CAPEm
Hc = sigC * alpc / Hm * sqrt(CAPEl)

P  = Hd + xis * Hs + xic*Hc

tem  = q + tstpi * (t1 + alpha2 * t2)
temteb = temteb_bar+teb - tem

D = m0 * ( 1 + mu * ( Hs -Hc )/qr01 ) * temteb
E = (tebs - teb  + sst_bar) / tau_e
Dryness = (temteb) / T0

# Stochastic Transition Rates

G7 = lambda x: 1- exp(x)

Rc= Matrix(zeros(4))

Rc[0,1] = G7(CAPEl) * G7(Dryness) / tau01
Rc[1,0] = G7(Dryness) / tau10
Rc[1,2] = G7(CAPE) * (1-G7(Dryness))/tau12
Rc[0,2] = G7(CAPE) * (1-G7(Dryness))/tau02
Rc[2,0] = (1 - G7(CAPE))  / tau20
Rc[2,3] = 1/tau23
Rc[3,0] = 1/tau30

for i in range(4):
    Rc[i,i] = - sum(Rc[i,:])

# Deterministic RHS

FT1 = P - Rad1
FT2 = Hc- Hs - Rad2
FTEB = (E - D) / Hb
FQ  = -tstpi  * P + D / Ht
FHS = alps * sigS * Hd / fdeq  -Hs

FDet = [ FT1, FT2, FTEB, FQ, FHS ]

# Stochastic RHS
sig0 = 1- sigC -sigD -sigS
sigs = Matrix([sig0, sigC, sigD, sigS])

FStoch = list(sigs.T * Rc)

FSRCE = FStoch[1:] + [FT1, FHS]

# Total RHS

FTot = FDet + FStoch

######################################################################
#                    Parameters RCE and linearization                #
######################################################################

param = dict(
    a1=0.5,
    a2=.5,
    a0=2,

    a1p=1.,  #Coefficient  of theta_eb in low CAPE, should be 1
    a2p=0.,  #Coefficient of q in low CAPE should be zero
    a0p=1.5,
    alpha2=.1, # gamma_2 in the papers (JASI,JASII,TAO,TCFD)
    alpha3=.1, # alpha_2 in the papers
    alpha4=2.,

    mu =2.0,

    tau01 = 1.0,
    tau02 = 3.0,
    tau10 = 1.0,
    tau12 = 1.0,
    tau20 = 3.0,
    tau30 = 5.0,
    tau23 = 3.0,

    #Units are K,
    theta_eb_m_theta_em = 11.0 / alpha_bar,
    theta_ebs_m_theta_eb = 10.0 / alpha_bar,
    moist0              = 30.0 / alpha_bar,
    thteb_st            = 0.0 / alpha_bar,


    # K/ day,
    qr01 = (1.0 / day_s) / ( alpha_bar / Te ),

    # J / kg
    CAPE0 = 400.0 / c**2,


    xis = 0.4,
    xic = 0.0,
    alpha_s = 0.25,
    alpha_c = 0.10

)

Hs, Hseq = symbols('Hs Hs_bar')
eps = symbols('epsilon')

rce_unknowns = [ m0, tau_e , capebar, qbar, qr02, Hseq, fdeq, fseq, fceq]

# d = solve([a.subs(param).subs(value_rce) for a in [ FT1, FT2, FTEB, FQ, FHS]], rce_unknowns)


phase_vars = [
        t1 ,
        t2 ,
        u1 ,
        u2,
        q,
        teb,
        sigC,
        sigD,
        sigS,
        Hs
        ]

linear_vars = [
        eps*t1 ,
        eps*t2 ,
        eps*u1 ,
        eps*u2,
        eps*q,
        eps*teb,
        fceq + eps*sigC,
        fdeq + eps*sigD,
        fseq + eps*sigS,
        Hseq + eps * Hs
        ]

RCE_params = [qbar, Hseq, fdeq, fseq, fceq, tau_e, capebar, qr02, m0]

for pat, sub in zip(phase_vars, linear_vars):

    FTot = [  x.subs(param).subs({pat:sub}).subs({eps:0}) for x in FTot]
    FStoch = [  x.subs(param).subs({pat:sub}).subs({eps:0}) for x in FStoch]
    FDet = [  x.subs(param).subs({pat:sub}).subs({eps:0}) for x in FDet]
    FSRCE = [  x.subs(param).subs({pat:sub}).subs({eps:0}) for x in FSRCE]




def linearize_expr(x, eps):
    lin = [xx.subs({eps:0}) for xx in x]
    o1  = [xx.series(eps, n=2).coeff(eps) for xx in x]

    return lin, o1
lin = [xx.subs({eps:0}) for xx in FTot]



# print "helo"
# FRCE, FPrime =  linearize_expr(FTot, eps)


class NMLFactory(object):
    dimensional = True

    tau01 = 1.0
    tau02 = 3.0
    tau10 = 1.0
    tau12 = 1.0
    tau20 = 3.0
    tau30 = 5.0
    tau23 = 3.0

    #Units are K
    theta_eb_m_theta_em = 11.0 / alpha_bar
    moist0              = 30.0 / alpha_bar


    # K/ day
    QR01 = (1.0 / day_s) / ( alpha_bar / Te )

    # J / kg
    CAPE0 = 400.0 / c**2


    xi_s = 0.4
    xi_c = 0.0
    alpha_s = 0.25
    alpha_c = 0.10

    # Heights

    # Rates
    Rc = np.zeros((4,4))


    # Outputs

    fceq = 0.0
    fseq = 0.0
    fdeq = 0.0

    capebar = 0.0


    def _calc_rates(self, CAPEb):
        c = self
        c.Rc = np.zeros((4,4))
        G7 = lambda x: 1.0 - np.exp(-x)

        D  = c.theta_eb_m_theta_em / c.moist0
        C  = CAPEb / c.CAPE0

        F = np.zeros(4)

        c.Rc[0,1] = G7(C) * G7(D) / c.tau01
        c.Rc[1,0] = G7(D) / c.tau10
        c.Rc[1,2] = G7(C) * (1.0-G7(D))/c.tau12
        c.Rc[0,2] = G7(C) * (1.0-G7(D))/c.tau02
        c.Rc[2,0] = (1.0 - G7(C))  / c.tau20

        c.Rc[2,3] = 1.0/c.tau23
        c.Rc[3,0] = 1.0/c.tau30

        for i in range(4):
            c.Rc[i,i] = - c.Rc[i,:].sum()

    def _rce_root(self, x):
        c = self

        sigC=x[ 0 ];
        sigD=x[ 1 ];
        sigS=x[ 2 ];
        CAPEb=x[ 3 ];

        c._calc_rates(CAPEb)

        sigs= np.zeros(4)
        sigs[1:] = x[:3]
        sigs[0]  = 1.0 -sigs.sum()

        F = np.zeros(4)

        F[:3]  = np.einsum('j,ji', sigs, c.Rc)[1:]

        Q=np.sqrt(CAPEb)/Hm;

        F[3]=c.QR01-sigD*Q-c.xi_s*c.alpha_s*sigS*Q-c.xi_c*c.alpha_c*sigC*Q;


        return F

    def _calc_rce_fracs(c):
        c.fceq, c.fdeq, c.fseq, c.capebar  = fsolve(c._rce_root, [.1, .1, .1, .1], factor=.1)


    def summary(c):
        print(c.fceq, c.fdeq, c.fseq, c.capebar)
        return np.array([c.fceq, c.fdeq, c.fseq, c.capebar])

    def print_namelist(c):
        raise NotImplementedError


def test_calc_rce_fracs():
    n = NMLFactory()
    n._calc_rce_fracs()
    out = n.summary()
    err=  out - np.array([0.00583515572402, 0.00209185536618, 0.0034864256103, 0.000960033827461])
    assert np.abs(err).sum() <  1e-5

    return 0


def test_rce_root():
    n = NMLFactory()
    out = n._rce_root([.1, .1, .1, .1])
    err = out -np.array([  0.30783428,   0.4680254 ,   0.11106667, -10.40866593])
    assert np.abs(err).sum() < 1e-5

    return 0

def test_sympy_matches_class():
    from numpy.testing import assert_almost_equal
    f = lambdify( rce_unknowns, FTot)
    ff = lambda x : f(*x)
    sympy_x = fsolve(ff, [.1]*9, factor=.1)
    print(sympy_x)

    n  = NMLFactory()
    n._calc_rce_fracs()
    class_x = [n.fceq, n.fdeq, n.fseq, n.capebar ]

    assert_almost_equal(sympy_x[:-1], class_x)






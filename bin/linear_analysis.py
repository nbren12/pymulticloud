
"""
Linear analysis


June 12, cross checked with KM08 code. The B MAtrices are the same

"""
import numpy as np
from numpy import exp, sqrt, pi, zeros
import matplotlib.pyplot as plt

hour = 3600
day  = 86400
km   = 1000

T = 2.9575e+04
L = 1.4827e+06
c = L / T
alpha_bar = 15.34


def grcpstruct(A, B, k_zonal):
    """
    Get growth rate, phase speed, and wave structure
    """
    scal = (2*np.pi * L / 4e7)
    k =  k_zonal * scal
    C =  - 1.0j * A * k + B(k)

    diag, V  = np.linalg.eig(C)

    gr = +np.real(diag)
    gr= gr / T*day

    if (abs(k - 0) > 1e-10):
        cp= - np.imag(diag)
        cp= cp / k * c
    else:
        cp = [ 0.0 ] * diag.shape[0]


    return gr, cp, V

def cp_growthrate_plots(A, B):
    scal = 2 * np.pi * L / 4e7 # Nondimensional wave number 1

    k_zonal = np.arange(0, 200 , 1)
    ks = k_zonal* scal         # List of wave numbers


    print("Plotting")

    Gr = []
    Ph = []

    for k in k_zonal:
        gr, cp, V = grcpstruct(A, B, k)

        Gr.append(gr)
        Ph.append(cp)

    Gr = np.vstack( Gr )
    Ph = np.vstack( Ph )

    print("Maximum growth rate occurs at k={} rate= 1/{:.2f}".format(Gr.max(axis=1).argmax(), 1/Gr.max()))

    xloc = plt.MaxNLocator(20)

    plt.figure(1)
    tmpx  = Ph * 0 + k_zonal[:,None]
    plt.scatter(tmpx, Ph, s=2.0)
    plt.grid('on')
    plt.gca().xaxis.set_major_locator(xloc)

    plt.figure(2)
    plt.scatter(tmpx, Gr, s=2.0)
    plt.gca().xaxis.set_major_locator(xloc)

    Ph = Ph[Gr > 0]
    tmpx = tmpx[Gr > 0]
    plt.grid('on')

    plt.figure(1)
    plt.scatter(tmpx, Ph, s= 30, c='r', edgecolor='w')

def x_deriv_mat(lmd_tld = 0.8, q_tld=0.9105):
    A = np.zeros((8,8))
    A[0,2]= -1.0 # u1
    A[1,3] = -1.0 #u2
    A[2,0] = -1.0 #t1
    A[3,1] = -1.0/4.0 #t2

    # Moisture equation
    A[4,0] = q_tld
    A[4,1] = q_tld * lmd_tld

    return A

def source_matrix():
    """
    Generate linearized matrix for MK08 source terms
    """
    alpha = .1
    tau_conv =  2 * hour / T
    tau_r    = 50 * day / T
    tau_c    =1 *hour /T
    tau_s    = 3 * hour /T

    a0 = 5
    a1 = .5
    a2 = .5
    a0p = 1.5  # This is an discrpeancy between code and paper
    gamma2 = .1
    gamma2p = 1.5
    alpha2 = 0.1

    lmbdst = 0.0

    thetap = 20 / alpha_bar
    thetam = 10 / alpha_bar

    theta_ebs_m_theta_eb = 10.0 / alpha_bar
    theta_eb_m_theta_em = 11.0 / alpha_bar

    xic = -0.5
    xis = 0.4

    # xis =.5
    # xic = 1.25



    alpha_s = .25
    alpha_c = 0.1

    mu = 0.35
    mu = 3
    hcdowndraft=0


    qr01 = 1.0 /day * T / alpha_bar
    ud =   0.1229


    def theta_em( t1, t2, q, alpha2 ):

        return q + 2 * sqrt( 2 ) / pi  * ( t1 + alpha2 * t2)

    def Qdp( t1, t2 ,teb, q,  tau_c0, a1, a2, a0, gam2):

        return ( a1 * teb + a2 * q - a0 * ( t1 + gam2 * t2 ) ) / tau_c0

    def Qdc( t1, t2, teb, a0p, gam2p, tau_conv ):

        return ( teb - a0p *( t1 + gam2p * t2 ) ) /tau_conv


    def Hdp( lmd, qdp, lmdbar, lmbdst, qbar):

        return ( ( 1.0 - lmdbar ) * qdp - qbar * lmd   )/ ( 1- lmbdst )


    zt = 15.75e3 / L
    zb = 500. / L

    u1, u2, t1, t2, teb, q, hc, hs = np.eye(8)



    # Deep heating

    Qd =  Qdp( t1, t2 ,teb, q,  tau_conv, a1, a2, a0, gamma2)
    Qc =  Qdc( t1, t2, teb, a0p, gamma2p, tau_conv )

    tem   = theta_em(t1, t2, q, alpha2)

    if theta_eb_m_theta_em > thetap:
        lmd_bar = 1.0
    elif theta_eb_m_theta_em < thetam:
        lmd_bar = lmbdst
    else:
        lmd_bar  = lmbdst + ( theta_eb_m_theta_em - thetam ) /  ( thetap - thetam)
    lmdp     = (teb - tem) /( thetap -thetam )

    qbar  = qr01 / ( ( 1- lmd_bar ) / (1 -lmbdst) * ( 1 + xis * alpha_s )              + xic *alpha_c * ( lmd_bar - lmbdst )/ ( 1 - lmbdst ) )


    Hd = Hdp( lmdp, Qd, lmd_bar, lmbdst, qbar)

    # Downdrafts

    hdbar = (1.0 - lmd_bar) / ( 1 - lmbdst ) * qbar
    hsbar = (1.0 - lmd_bar) / ( 1 - lmbdst ) * qbar * alpha_s
    hcbar = (lmd_bar - lmbdst) / ( 1 - lmbdst ) * qbar * alpha_c

    dbar  = qr01 * 2 * sqrt( 2 ) / pi * zt

    m0    = dbar / ( 1 + mu / qbar * ( hsbar - hcdowndraft * hcbar ) ) / (theta_eb_m_theta_em)

    D =  m0 * (mu / qbar *( hs - hcdowndraft*hc ) ) * theta_eb_m_theta_em         + m0 * ( 1 + mu/qbar * ( hsbar - hcbar ) ) * ( teb - tem )

    # Evaporation
    tau_e = theta_ebs_m_theta_eb  * zb / dbar

    E = -teb / tau_e

    # Total precip

    Pr = Hd + xis * hs + xic * hc



    B = np.zeros(( 8, 8 ))

    B[0,:] =  - u1 * ud
    B[1,:] =  - u2 * ud
    B[2,:] = Pr - t1/tau_r
    B[3,:] = hc - hs - t2 / tau_r
    B[4,:] = E - D / zb
    B[5,:] = -2 * sqrt(2) /pi * Pr + D / zt
    B[6,:] = ( alpha_c / ( 1.0 - lmbdst) *  ( lmdp * qbar + ( lmd_bar -lmbdst  )* Qc ) - hc  ) / tau_c
    B[7,:] = ( alpha_s * Hd - hs  ) / tau_s

    Bfun = lambda k : B


    return Bfun

A  = x_deriv_mat()
B  = source_matrix()


print("hello")
cp_growthrate_plots(A, B)
plt.show()

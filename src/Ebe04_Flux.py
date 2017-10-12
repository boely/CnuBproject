from math import *
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.integrate as integrate

# ============================================================================
# Ebe04: 
#   Eberle 2004
#   Relic Neutrino Absorption Spectroscopy
#   arxiv:hep-ph/0401203
# ============================================================================

# =================
# "Basic formula's"
#  - Resonance energy
#  - Hubble parameter
# =================

def E_resonance(m_n):
    m_z = 91.2e9                    # [eV]
    return m_z**2 / (2. * m_n)                                                  # (eq 1)

def Hubble(z, h=0.71, Om=0.3, Ol=0.7, Ok=0.):
    H0 = 100 * h                    # [km/s/Mpc]
    return sqrt(H0**2 * ( Om * (1+z)**3 + Ok * (1+z)**2 + Ol) )                 # (eq 8)


# =================
# Survival Probability
# =================

# -- as a function of E/Eres and z
def Ebe04_survival_probability_P(E_frac_Eres, z, 
                                 h=0.71, Om=0.3, Ol=0.7, Ok=0.):
    '''survival probability function
    returns survival probability for the E_frac_Eres values within range:
    1/(1+z) < E_frac_Eres < 1 ''' 

    ann_prob = 0.71/h * 0.03        # anihilation probability                   # (eq16)
    assert(float(Om + Ol + Ok) == 1.), "Not My Universe"

    if E_frac_Eres < 1./(1+z) or E_frac_Eres > 1.:
        return 1.                   # No absorption
    else:
        return exp(-ann_prob * ((1/E_frac_Eres**3) / sqrt(Om*(1/E_frac_Eres**3) 
        + Ok*(1/E_frac_Eres**2) + Ol)))                                         # (eq15)


# =================
# neutrino flux
# =================

def L_remained(z, z_max, n_min_alpha):                                          # what is remained of emisivity after deviding out factors independent of z
    z_min = 0
    if z > z_min:
        if z_max > z:
            return (1 + z) ** n_min_alpha
        else: return 0.
    return 0.


def Ebe04_neutrino_flux_earth_F(E_frac_Eres, z_max, n_min_alpha):
    '''Flux per neutrino flavor '''

    # - function to integrate for Flux
    f1_to_integrate = lambda z: 1/Hubble(z) * Ebe04_survival_probability_P(E_frac_Eres, z) * L_remained(z, z_max, n_min_alpha)
    integrant1, err1 = integrate.quad(f1_to_integrate, 0, np.inf)

    # - function to integrate for Flux in case of no absorption (P == 1)
    f2_to_integrate = lambda z: 1/Hubble(z) * L_remained(z, z_max, n_min_alpha)
    integrant2, err2 = integrate.quad(f2_to_integrate, 0, np.inf)

    if integrant2 != 0:
        # return fraction F/F_no abs
        return integrant1 / integrant2
    else:
        return 1


# =================
# Plot functions
# =================

def plot_P():
    redshift = [2., 5., 20.]
    for z in redshift:
        E_frac_Eres = np.linspace(1/(1-z), 1.5, 1000)
        P = []
        for i in range(len(E_frac_Eres)):
            P.append(Ebe04_survival_probability_P(E_frac_Eres[i], z))
        labeltxt = "z = %i" %z
        plt.plot(E_frac_Eres, P, "-", ms = 5, label = labeltxt)
    plt.legend(loc = 'lower right')
    plt.xscale('log')
    plt.xlabel("E/Eres")
    plt.ylabel("P(E(1+z),z)")
    plt.xlim(0.01, 1.5)
    plt.ylim(0, 1.1)
    plt.title('Survival Probability, Ebe04')
    plt.savefig("Survival_Probability_Ebe04.png")
    plt.show()

def plot_F():
    redshifts_max = [2, 5, 20]
    n_min_alphas = [0, 2, 4]

    for z_max in redshifts_max:
        if z_max == 2: col = 'r'
        elif z_max == 5: col = 'b'
        elif z_max == 20: col = 'g'

        for n_min_alpha in n_min_alphas:
            if n_min_alpha == 4: linst = '-'
            elif n_min_alpha == 2: linst = '--'
            elif n_min_alpha == 0: linst = ':'

            E_frac_Eres = np.logspace(-2, 0, 1000)  # from 10^-2 to 10^0 in 1000 logarithmic steps
            F_frac_Fnoabs = []
            for e in E_frac_Eres:
                F_frac_Fnoabs.append(Ebe04_neutrino_flux_earth_F(e, z_max, n_min_alpha))
            labeltxt = "z = %i n-alpha = %i" %(z_max, n_min_alpha)
            print 'plot', labeltxt
            plt.plot(E_frac_Eres, F_frac_Fnoabs, linestyle = linst, color = col, label = labeltxt)

    plt.legend(loc='best')
    plt.title("Flux, Ebe04")
    plt.xlim(1e-2, 1)
    plt.xlabel("E/E_res")
    plt.xscale('log')
    plt.ylim(0, 1)
    plt.ylabel("F/F_no_abs")
    plt.savefig("Ebe04_Flux.png")
    plt.show()

# =================
# main
# =================

def main():
    # -- plot absorption probability as a function of E/Eres
    plot_P()

    # -- plot Flux as a fraction of Flux in case of no absorbtion, as a function of E/Eres
    plot_F()


if __name__ == "__main__":
    sys.exit(main())
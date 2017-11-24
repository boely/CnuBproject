from math import *
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.integrate as integrate
import os
import ConfigParser

# =================
# Plot functions
# =================
sys.path.append("src/.")
from NeutrinoFlux import *

def plot_sigma(flux):
    flux.m_n = 0.1
    flux.E_resonance()
    E = np.logspace(-2, 3, 10000)*flux.Eres
    
    plt.plot(E, np.array([Oli05_crossection_Sigma(x) for x in E]), '-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel("E$\\nu$ [eV]", size =14, position=(1, 0), ha='right')
    plt.ylabel("$\sigma$ [nb]", size =14, position=(0,1), ha='right')
    plt.show()

def plot_P(flux):
    redshift = [2., 5., 20.]
    flux.m_n = 0.1
    flux.E_resonance()
    for z in redshift:
        if z == 2: col = 'r'
        elif z == 5: col = 'b'
        elif z == 20: col = 'g'
        flux.z_max = z
        E = np.linspace(flux.Eres/(1-z), flux.Eres*1.5, 1000)
        P = []
        for e in E:
            P.append(flux.Ebe04_survival_probability_P(e))
        labeltxt = "z = %i" %z
        plt.plot(E, P, "-", ms = 5, color = col, label = labeltxt)
    plt.legend(loc = 'lower right')
    plt.xscale('log')
    plt.xlabel("E (eV)")
    plt.ylabel("P(E(1+z),z)")
    plt.xlim(0.01*flux.Eres, 1.5*flux.Eres)
    plt.ylim(0, 1.1)
    plt.title('Survival Probability, Ebe04')
    plt.text(2e22, 0.25, r'$ m_{\nu} = 0.1 eV $', size = 18)
    plt.show()


def plot_F(flux, frac = False, Z_decay = False):
    '''Function to plot the flux.
    If frac = True, the flux will be plotted divided by the flux when there
    is no absorption.'''
    # Fod03, fitting procedure gives the most probable values for Emax, alpha and n.
    redshifts_max = [20]  #[2, 5, 20]                                           # z_max
    alphas        = [2] #[0.5, 1, 1.5, 2]                                       # spectral index -> (Ebe04: [1-2])
    ns            = [6] #[2, 4, 6]                                              # powers of (1+z) in activity Ebe04 eq27

    for z_max in redshifts_max:
        if z_max == 2: col = 'r'
        elif z_max == 5: col = 'b'
        elif z_max == 20: col = 'g'

        for n in ns:
            for alpha in alphas:
                if n - alpha == 4: linst = '-'
                elif n - alpha == 2: linst = '--'
                elif n - alpha == 0: linst = ':'
                else: continue # <= for now only analyse (n-alpha)'s from fig4

                flux.z_max = z_max
                flux.n = n
                flux.alpha = alpha
                flux.m_n = 0.1
                flux.E_resonance()
                E = np.logspace(-2, 3, 10000)*flux.Eres           # from 10^a to 10^b in c logarithmic steps
                F = []
                F1 = []

                if frac == True:

                    Fnoabs = []
                    F_frac_Fnoabs = []
                    F1noabs = []
                    F1_frac_Fnoabs = []

                    for e in E:
                        F.append(flux.Ebe05_neutrino_flux_earth_F(e, Z_decay, absorption = True))
                        Fnoabs.append(flux.Ebe05_neutrino_flux_earth_F(e, Z_decay, absorption = False))
                        F_frac_Fnoabs.append(F[-1]/Fnoabs[-1])

                        if Z_decay == True:
                        # If flux + secondary flux is plotted: Z_decay = True,
                        # add plot of primary flux only in black, to compare:
                            F1.append(flux.Ebe05_neutrino_flux_earth_F(e, Z_decay=False, absorption = True))
                            F1noabs.append(flux.Ebe05_neutrino_flux_earth_F(e, Z_decay=False, absorption = False))
                            F1_frac_Fnoabs.append(F1[-1]/F1noabs[-1])

                    labeltxt = "z_max = %i, n = %i, alpha = %i" %(z_max, n, 
                                                                  alpha)
                    print 'plot', labeltxt
                    plt.plot(E/flux.Eres, F_frac_Fnoabs, linestyle = linst, 
                             color = col, label = labeltxt)
                    if Z_decay == True: 
                        plt.plot(E/flux.Eres, F1_frac_Fnoabs, linestyle = linst, 
                             color = 'black', label = labeltxt)
                    plt.ylabel("F/F_no_abs")
                    plt.text(5e20, 0.65, r'$ m_{\nu} = 0.1 eV $', size = 18)
                    #plt.ylim(0, 1.1)
                    plt.xlim(1e-2, 1)
                    plt.xlabel("E/Eres")


                else:

                    for e in E:
                        F.append(flux.Ebe05_neutrino_flux_earth_F(e, Z_decay, absorption = True) * e**2)     # < ==  F * E^2

                    labeltxt = "z_max = %i, n = %i, alpha = %i" %(z_max, n, 
                                                                  alpha)
                    print 'plot', labeltxt
                    plt.plot(E, F, linestyle = linst, color = col, 
                             label = labeltxt)
                    plt.ylabel("F * E^2")
                    plt.yscale('log')
                    plt.text(5e21, 1.7e-46, r'$ m_{\nu} = 0.1 eV $', size = 18)
                    #plt.xlim(1e-2 * Eres, Eres)
                    plt.xlabel("E (eV)")

    plt.legend(loc = 'best')
    if Z_decay == True:
        plt.title("Flux + Z-decay flux, Ebe04+Ebe05")
    elif Z_decay == False:
        plt.title("Flux, Ebe04")
    plt.xscale('log')
    plt.show()



def main():

    args = sys.argv[1:]
    if "-h" in args or "--help" in args or len(args) < 2:
        usage()
        sys.exit(2)

    config = ConfigParser.RawConfigParser()
    config.read(args[0])
    sys.stderr.write('Reading config file %s \n'% args[0])

    flux = NeutrinoFlux()
    flux.z_max = config.getfloat("Flux", "z_max")
    flux.alpha = config.getfloat("Flux", "alpha") 
    flux.n     = config.getfloat("Flux", "n") 
    flux.m_n   = config.getfloat("Flux", "m_n") 
    flux.eta0  = config.getfloat("Flux", "eta0") 
    flux.j     = config.getfloat("Flux", "j") 
    flux.Eresolution = config.getfloat("Telescope", "energy_resolution")

    prefix_outfilename = args[1]
    flux = NeutrinoFlux()
    
    # -- plot crossections
#    plot_sigma()

    # -- plot absorption probability as a function of Energy
#    plot_P()

    # -- plot Flux as a fraction of Flux in case of no absorbtion, 
    #    as a function of Energy/Eres
    plot_F(flux, frac = True, Z_decay = False)
    plot_F(flux, frac = True, Z_decay = True)




if __name__ == "__main__":
    sys.exit(main())


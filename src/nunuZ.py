import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime as dt
import math
from scipy.constants import *

# globals
GeV = 1E9
eV = 1

class nunuZ:
    def __init__(self, neutrino_mass = 1e-5*eV):
        self.neutrino_mass = neutrino_mass
        self.sinthetaW2 = 0.2223
        self.Gamma_Z = 2.4952 * GeV
        self.Massa_Z = 91.1876 * GeV
        self.G_f = physical_constants["Fermi coupling constant"][0] # in eV^-2
        self.G_f2 = self.G_f**2 
        
        self.mbarn_per_GeV2 = 0.389379
        self.mbarn_per_eV2 = 0.389379e-18
        self.mbarn_to_cm2 = 1e-27
        self.final_state_particle = None
        
    def F2(self, Y):
        assert Y !=0
        return (3*Y**2 + 2*Y - 2*(1+Y)**2 * math.log(1+Y) )/Y**3

    def F1(self, Y):
        assert Y != 0
        return ( Y**2 + 2*Y -2*(1+Y)*math.log(1+Y))/Y**3

    def set_neutrino_mass(self, neutrino_mass):
        self.neutrino_mass = neutrino_mass
    
    def get_cms_energy(self):
        assert self.neutrino_mass > 0
        return (0.5*self.Massa_Z**2/self.neutrino_mass, "eV")
        
    def set_final_state_particle(self, particle):
        self.final_state_particle = particle

    def nunub_ffb_s(self, energy):
        # equation quigg, 19  p6
        assert self.final_state_particle
        (Nc, Qf, tau3) = self.final_state_particle
    
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2
        Lf = tau3 - 2*Qf * self.sinthetaW2
        Rf = -2*Qf * self.sinthetaW2
        LR_factor = (Lf**2 + Rf**2)/( (1-Y)**2 + (self.Gamma_Z/self.Massa_Z)**2)

        # switch to GeVs
        return self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) * Nc/(3*pi) * \
            LR_factor * self.mbarn_per_GeV2 * self.mbarn_to_cm2

    def nunub_nunub_t(self, energy):
        # equation quigg, 20  p7
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2

        # switch to GeV
        return (self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / pi ) * self.F1(Y) *\
            self.mbarn_per_GeV2 * self.mbarn_to_cm2

    def nunub_nunub_st(self, energy):
        # equation quigg, 21  p7
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2
        factor = (Y-1) / ( (1-Y)**2 + (self.Gamma_Z/self.Massa_Z)**2)

        # switch to GeV, Y is dimensionless
        return (self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / (2*pi) ) * self.F2(Y) *\
            factor * \
            self.mbarn_per_GeV2 * self.mbarn_to_cm2

    def nunu_nunu_t(self, energy):
        # equation quigg, 24  p7
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2

        # switch to GeV, Y is dimensionless
        return (self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / pi ) * (1/(1+Y)) * \
            self.mbarn_per_GeV2 * self.mbarn_to_cm2

    def nunu_nunu_u(self, energy):
        # equation quigg, 25  p7
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2

        # switch to GeV, Y is dimensionless
        return (self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / pi ) * \
            ( (1/(1+Y)) + math.log(1+Y) /( Y*(1+0.5*Y)) ) * \
            self.mbarn_per_GeV2 * self.mbarn_to_cm2

    def nunub_llb_t(self, energy):
        # equation quigg, 22  p7
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2

        # switch to GeV
        return (4*self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / pi ) * self.F1(Y) * \
            self.mbarn_per_GeV2 * self.mbarn_to_cm2

    def nunub_llb_st(self, energy):
        # equation quigg, 23  p7
        Y = 2 *self.neutrino_mass * energy / self.Massa_Z**2
        factor = (Y-1) / ( (1-Y)**2 + (self.Gamma_Z/self.Massa_Z)**2)

        # switch to GeV, Y is dimensionless
        return (4*self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / (2*pi) ) * self.F2(Y) * \
            (self.sinthetaW2 - 0.5) * factor * \
            self.mbarn_per_GeV2 * self.mbarn_to_cm2

    
    def nunub_ffb_s_channel(self, UHE_energy, particle):
        self.set_final_state_particle(particle)
        return self.nunub_ffb_s(UHE_energy)
                          
    def nunub_nunub_t_channel(self, UHE_energy):
        return self.nunub_nunub_t(UHE_energy)

    def nunub_nunub_st_channel(self, UHE_energy):
        return self.nunub_nunub_t(UHE_energy)

    def nunu_nunu_t_channel(self, UHE_energy):
        return self.nunu_nunu_t(UHE_energy)

    def nunu_nunu_u_channel(self, UHE_energy):
        return self.nunu_nunu_u(UHE_energy)

    def nunub_llb_t_channel(self, UHE_energy):
        return self.nunub_llb_t(UHE_energy)

    def nunub_llb_st_channel(self, UHE_energy):
        return self.nunub_llb_st(UHE_energy)
    
def main():
    
    neutrino_mass = 1e-5
    nnZ = nunuZ(neutrino_mass)
    print nnZ.get_cms_energy()
    fig = plt.figure()

    energy = np.arange(1E25, 1E28, 1E24) #eV

    # particle: Nc, charge, isospin
    quark = (3, 1/3., 1.)
    neutrino = (1, 0, 1.)
    
    Xsection = []
    Xsection.append( np.asarray( [nnZ.nunub_ffb_s_channel(e, neutrino) for e in energy]) )
    Xsection.append( np.asarray( [nnZ.nunub_ffb_s_channel(e, quark) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunub_nunub_t_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunub_nunub_st_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunu_nunu_t_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunu_nunu_u_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunub_llb_t_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunub_llb_st_channel(e) for e in energy] ) )


    Xsec = Xsection[0]
    for i in range(1,len(Xsection)):
        Xsec += Xsection[i]

    ax1 = fig.add_subplot(111)
    
    for i, X in enumerate(Xsection):
        ax1.loglog(energy, X)
    ax1.loglog(energy, Xsec, label="total")
    ax1.set_ylim([5E-36, 1E-30])
    ax1.set_xlim([1E25, 1E28])
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

import numpy as np
from math import *

# natures constants
from scipy.constants import codata
# particle data
try:
    import pypdt
    tbl = pypdt.ParticleDataTable()
    GZ = tbl[23].width # in GeV
    Mz  = tbl[23].mass  # in GeV
except ImportError:
    print('\n Warning: Pypdt module not installed')
    GZ = 2.4952
    Mz = 91.1876

class cross_section():
    def __init__(self):
        self.GF = codata.value('Fermi coupling constant') # in GeV^-2
        self.GammaZ = GZ
        self.MassZ = Mz
        self.p = 0.2e-9 # in GeV
        self.m = 0.1e-9 # in GeV
    
    def sigma_r(self, Eprime, s):
        # Lunardini et al (2013), equation A.4
        xi = self.GammaZ**2/self.MassZ**2
        return ( (self.GF * self.GammaZ * self.MassZ) / \
            (sqrt(2.) * Eprime * sqrt(self.p**2 + self.m**2)) * \
            ( (s/(1 + xi)) \
              - (self.MassZ**2 * (xi -1) / (sqrt(xi) * (1 + xi)**2)) * \
              atan( ((1. + xi) * s - self.MassZ**2) / (self.MassZ**2 * sqrt(xi)))
               + \
               (self.MassZ**2 / (1 + xi)**2) * log( (1+xi) * s**2 \
                                                    - 2.*self.MassZ**2 \
                                                    + s*self.MassZ**4) \
        ))

             
    def sigma_resonant(self, Eprime):
        # Eprime = energy of the UHEnu
        s_plus = 2*Eprime * (sqrt(self.p**2 + self.m**2) + self.p)
        s_minus = 2*Eprime * (sqrt(self.p**2 + self.m**2) - self.p)
        return self.sigma_r(Eprime, s_plus) - self.sigma_r(Eprime, s_minus)


    def Sigma_r(self, K):
        # D'Olivo et al 2005, equation 22
        P = self.p
        Ep = sqrt(self.p**2 + self.m**2)
        xi = self.GammaZ**2/self.MassZ**2
        teller = (1+xi)*4*K**2*(Ep+P)**2 - 4*self.MassZ**2*K*(Ep+P) + self.MassZ**4
        noemer = (1+xi)*4*K**2*(Ep-P)**2 - 4*self.MassZ**2*K*(Ep-P) + self.MassZ**4

        aid1 = (2*K*(1+xi)*(Ep+P)-self.MassZ**2)/(self.GammaZ*self.MassZ)
        aid2 = (2*K*(1+xi)*(Ep-P)-self.MassZ**2)/(self.GammaZ*self.MassZ)
        

        # convert cross section units:
        # 1 mb = 1/2.56819 GeV-2, 1mb = 1e-27cm^2
        C = (1/2.56819)*1e-27
        
        return ( C * (2*sqrt(2.)*self.GF * self.GammaZ * self.MassZ) / (K * Ep) * \
                 (
                     (1/(1 + xi))
                     + self.MassZ**2/(4*K*P*(1+xi)**2) * log(teller/noemer) \
                     + ((1-xi)*self.MassZ**3 / ((1+xi)**2 * 4*K*P*self.GammaZ)) * ( atan(aid1) - atan(aid2))
                )
        )


    
import matplotlib.pyplot as plt

def main(argv):
    X = cross_section()
    X.p = 6.08e-4*1e-9
    X.m = 1e-4 * 1e-9

    plt.figure(99)
    E = np.arange(14,18,.0001)
    
    plt.semilogy(E, np.array([X.Sigma_r(pow(10,x)) for x in E]), '-')
    plt.xlabel("$\log$(E$\\nu$) [GeV]", size =14, position=(1, 0), ha='right')
    plt.ylabel("$\sigma$ [cm$^2$]", size =14, position=(0,1), ha='right')
    plt.show()

import sys
    
if __name__ == "__main__":
    main(sys.argv)

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
        self.sinthetaW2 = 0.231   
        self.Gamma_Z = 2.4952 * GeV
        self.Massa_Z = 91.1876 * GeV
        self.G_f = physical_constants["Fermi coupling constant"][0] # in eV^-2
        self.G_f2 = self.G_f**2 
        
        self.mbarn_per_GeV2 = 0.389379
        self.mbarn_per_eV2 = 0.389379e-18
        self.mbarn_to_cm2 = 1e-27
        self.final_state_particle = None
        self.verbose = True
        
        
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

    #---------------------------------------------
    #--/ individual cross-sections (Quigg, page 6/7)
    #---------------------------------------------
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
        return (4*self.G_f2 * (self.neutrino_mass/GeV) * (energy/GeV) / (pi) ) * self.F2(Y) * \
            (self.sinthetaW2 - 0.5) * factor * \
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


    #-------------------------------------------------------------------------------------
    #-- list of all processes (set final state particle type and call cross-section computation) 
    #-------------------------------------------------------------------------------------    
    def nunub_ffb_s_channel(self, UHE_energy, particle):
        self.set_final_state_particle(particle)
        return self.nunub_ffb_s(UHE_energy)
                          
    def nunub_nunub_t_channel(self, UHE_energy):
        return self.nunub_nunub_t(UHE_energy)

    def nunub_nunub_st_channel(self, UHE_energy):
        return self.nunub_nunub_st(UHE_energy)

    def nunub_llb_t_channel(self, UHE_energy):
        return self.nunub_llb_t(UHE_energy)

    def nunub_llb_st_channel(self, UHE_energy):
        return self.nunub_llb_st(UHE_energy)

    def nunu_nunu_t_channel(self, UHE_energy):
        return self.nunu_nunu_t(UHE_energy)

    def nunu_nunu_u_channel(self, UHE_energy):
        return self.nunu_nunu_u(UHE_energy)
            
    #-------------------------------------------------------------
    #-- print checks & details of cross-section computation & plotting
    #-------------------------------------------------------------

    # details of cross-section computation
    def printdetails(self, Xsection, energy):
       
          # find resonance  energy (and index)
          index_res = Xsection[0].argmax(max(Xsection[0]))
          print(" Resonance energy at %5.2g eV" % energy[index_res])

          # test (Z decay width)
          Xsec_max_s_channel = Xsection[0][index_res] + Xsection[1][index_res]  + Xsection[2][index_res] + Xsection[3][index_res]
          print("\n Test: s-channel breakdown (Z decay):")
          print("    Z decay width (nn) = %5.2f %%" % (100*Xsection[0][index_res]/Xsec_max_s_channel))
          print("    Z decay width (ll) = %5.2f %%" % (100*Xsection[1][index_res]/Xsec_max_s_channel))
          print("    Z decay width (qq) = %5.2f %%" % (100*(Xsection[2][index_res]+Xsection[3][index_res])/Xsec_max_s_channel))

          # print cross-section at the resonance
          print("\n Cross-section at resonance:")
          print("    Max (s-Z (qq))  = %6.3g" % (Xsection[2][index_res]+Xsection[3][index_res]))
          print("    Max (s-Z (nn))  = %6.3g" % Xsection[0][index_res])
          print("    Max (s-Z (ll))  = %6.3g" % Xsection[1][index_res])
          print ("    ---")
          print("    Max (t-Z)       = %6.3g" % Xsection[4][index_res])
          print("    Max (st-Z)      = %6.3g" % Xsection[5][index_res])
          print ("    ---")
          print("    Max (t-W-ll)    = %6.3g" % Xsection[6][index_res])
          print("    Max (st-ll)     = %6.3g" % Xsection[7][index_res])
          print ("    ---")
          print("    Max (t-nnnn)    = %6.3g" % Xsection[8][index_res])
          print("    Max (st-nnnn)   = %6.3g" % Xsection[9][index_res])

          return

    #-- produce plot of cross section versus energy (style: Quigg, fig 7
    def produceplot(self, Xsection, Xsec, energy):
 
           fig = plt.figure()
           plt.cla()
           ax1 = fig.add_subplot(111)
    
           # --total cross-section    
           ax1.loglog(energy, Xsec, label="total", color = 'black', linewidth = 2)
          # s-channel Z (19)
           ax1.loglog(energy, Xsection[0], label="neutrino", color = 'black', linewidth = 1)
           ax1.loglog(energy, Xsection[1], label="lepton", color = 'red', linewidth = 1)
           ax1.loglog(energy, Xsection[2]+Xsection[3],label="quarks", color = 'lightgreen', linewidth = 1)
           #-- t-channel and st-interference Z [nunubar final state] (20 and 21)
           ax1.loglog(energy, Xsection[4], label="t", color = 'blue', linewidth = 1)
           #####ax1.loglog(energy, Xsection[5], label="st", color = 'red', linewidth = 1)   # s-t channel interference (21)
           # -- t-channel and st-interference W  [ll final state] (22 and 23)
           ax1.loglog(energy, Xsection[6], label="t", color = 'orange', linewidth = 1)
           #####ax1.plot(energy, Xsection[7], label="st", color = 'red', linewidth = 1)     #s-t channel interference (23)
           #-- non-resonant (24 and 25)
           ax1.loglog(energy, Xsection[8], label="t", color = 'pink', linewidth = 1)
           ax1.loglog(energy, Xsection[9], label="t", color = 'lightblue', linewidth = 1)

           ax1.set_ylim([5E-36, 1E-30])
           ax1.set_xlim([1E26, 1E27])
           plt.grid(True)
           #plt.legend()
           plt.show()

           return



#=======
#== MAIN
#=======
def main():
    
    neutrino_mass = 1e-5
    nnZ = nunuZ(neutrino_mass)
    print("Resonance energy at %6.2g" % nnZ.get_cms_energy()[0])
    
    energy = np.arange(1E25, 1E28, 1E24) #eV

    #-- define particle types
    up_quark = (3, 2/3., 1)          
    down_quark = (3, -1/3., -1)   
    neutrino = (1, 0, 1)               
    lepton = (1, -1, -1)              
    
   #------------------------------------ 
   # [1] compute individual cross sections
   #------------------------------------ 
    Xsection = []
    #-- s-channel Z (19)
    Xsection.append( np.asarray( [nnZ.nunub_ffb_s_channel(e, neutrino) for e in energy]) )
    Xsection.append( np.asarray( [nnZ.nunub_ffb_s_channel(e, lepton) for e in energy]) )
    Xsection.append( np.asarray( [nnZ.nunub_ffb_s_channel(e, up_quark) for e in energy] ) )     
    Xsection.append( np.asarray( [nnZ.nunub_ffb_s_channel(e, down_quark) for e in energy] ) ) 
   #-- t-channel and st-interference Z [nunubar final state] (20 and 21)
    Xsection.append( np.asarray( [nnZ.nunub_nunub_t_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunub_nunub_st_channel(e) for e in energy] ) )
     #-- t-channel and st-interference W  [ll final state] (22 and 23)
    Xsection.append( np.asarray( [nnZ.nunub_llb_t_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunub_llb_st_channel(e) for e in energy] ) )
    #-- neutrino-neutrino (24 and 25)
    Xsection.append( np.asarray( [nnZ.nunu_nunu_t_channel(e) for e in energy] ) )
    Xsection.append( np.asarray( [nnZ.nunu_nunu_u_channel(e) for e in energy] ) )

    #----------------------------------------------------
    # [2] compute final cross-section (s-channel summation)
    #----------------------------------------------------
    # implement number of families for s-channel decay products
    Nfamilies = 3
    for i in range(0,len(Xsection[0])):    
         Xsection[0][i] = Nfamilies*Xsection[0][i]         # neutrino's
         Xsection[1][i] = Nfamilies*Xsection[1][i]         # leptons
         Xsection[2][i] = Nfamilies*Xsection[2][i]         # up-quarks
         Xsection[3][i] = (Nfamilies-1)*Xsection[3][i]  # down-quarks (ttbar excluded)

    Xsec = Xsection[0].copy()
    for i in range(1,len(Xsection)):    
        Xsec += Xsection[i]

    #------------------------------
    # [3] print details and show plot
    #------------------------------
    if nnZ.verbose:
        nnZ.printdetails(Xsection,energy)
        nnZ.produceplot(Xsection, Xsec, energy)





if __name__ == "__main__":
    main()

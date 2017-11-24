# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:59:51 2017

@author: mullerrs
"""
import os
import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from scipy.interpolate import interp1d
import scipy.integrate as integrate



# =================
# number of events
# =================

def events(Flux):
    '''Function to convert flux to amount of events'''
    Detector_eff = 0.01 # [fraction] for now assume 1% detector efficiency       <== Gauss
    t = 10              # [years]
    A = pi*500**2       # [m^2]      pi * (500m)^2 = opp 1 Antares blok
    Theta = 2*pi        # [sr]       Full sphere is 4pi sr
                        #            --> here half of the angular resolution
    return (60*60*24*365.25 * t) * A * Theta * Detector_eff * Flux


# =================
# Flux data
# =================

def spectrum_func(e, H, E_H0, FE2_H0, E_H1, FE2_H1 ):    # mean expected events
    '''Function that returns a function for the null hypothesis H = H0 
    (constant spectrum), and H = H1 hypothesis (dip in spectrum) '''

    # interpolate the separate datapoints
    f_H0 = interp1d(E_H0, FE2_H0)
    f_H1 = interp1d(E_H1, FE2_H1)

    if H == 'H0':
        return f_H0(e)
    elif H == 'H1':
        return f_H1(e)


def fill_hist( func, H, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Function to plot a histogram according to a function func, given the 
    hypothesis H'''
    bins = 100
    h1 = ROOT.TH1D( 'h1','hypotheses', bins, min(E_H0), max(E_H0))

    Binwidth = (max(E_H0)-min(E_H0))/bins

    for be in range( 1, h1.GetNbinsX() + 1 ):
        e = h1.GetXaxis().GetBinCenter( be )
        b = h1.GetBin( be ) # global bin number
        z = func(e, H, E_H0, FE2_H0, E_H1, FE2_H1 )
        h1.SetBinContent( b, z )

    return h1


def draw_H0_H1(E_H0, FE2_H0, E_H1, FE2_H1):
    '''Function that draws H0 and H1 from imported data in a histogram'''
    c1 = ROOT.TCanvas()
    h1 = fill_hist( spectrum_func , 'H0', E_H0, FE2_H0, E_H1, FE2_H1)
    h1.SetLineStyle(3)
    h1.Draw('LSAME')
    h2 = fill_hist( spectrum_func , 'H1', E_H0, FE2_H0, E_H1, FE2_H1)
    h2.Draw('LSAME')
    h2.SetXTitle("Energy (eV)")
    h2.SetYTitle("Flux per bin")
    h2.SetTitle("Search for dip")
    c1.Update()
    c1.SaveAs("original_functions_H0_and_H1.png")
    h1.Delete()
    h2.Delete()


# =================
# Accept-reject function
# =================

def accept_reject(E_H0, FE2_H0, E_H1, FE2_H1, H ):
    '''Function for acceptence, rejection in the ranges E_min-E_max, 0-FE2_max
    an accepted energy value will be returned'''
    if H == 'H1': 
        E = E_H1
        FE2 = FE2_H1
    elif H == 'H0':
        E = E_H0
        FE2 = FE2_H0

    # bounds of the box to create random numbers in:
    E_min = min(E)
    E_max = max(E)
    FE2_max = max(FE2)

    rand_e_xi = ROOT.gRandom.Rndm() * (E_max - E_min) + E_min                   # uniformly distributed number between E_min and E_max
    rand_f_ui = ROOT.gRandom.Rndm() * FE2_max                                   # uniformly distributed number between 0 and FE2_max

    #  -  accept-reject function
    while(spectrum_func( rand_e_xi, H, E_H0, FE2_H0, E_H1, FE2_H1) <= rand_f_ui ):
        rand_e_xi = ROOT.gRandom.Rndm() * (E_max - E_min) + E_min
        rand_f_ui = ROOT.gRandom.Rndm() * FE2_max
    return rand_e_xi


# =================
# Likelihood
# =================

def log_likelihood(hist, H, N_events, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Function that returns the log likelihood for a given hypothesis 'H' with
    respect to the given histogram 'hist' with your data. 
    NB: make sure FE2_H0 and FE_H1, interpolated are normalized to 1'''
    loglik = 0.

    # -- loop over bins
    for i_bin in range( 1, hist.GetNbinsX() + 1 ):
        e = hist.GetBinCenter( i_bin )                                          # energy (centre of bin) 
        mu_bin = spectrum_func(e, H, E_H0, FE2_H0, E_H1, FE2_H1) * hist.GetBinWidth(1) * N_events     # theoretical amount of counts in bin
        Nevt_bin = hist.GetBinContent( i_bin )                                  # the amount of counts in bin found in the data of histogram h

        if ROOT.TMath.Poisson( Nevt_bin, mu_bin ) > 0:                          # <= check if&when this will happen
            loglik += ROOT.TMath.Log( ROOT.TMath.Poisson( Nevt_bin, mu_bin ) )

    return loglik


def LLR(hist, N_events, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Function that caluclates the likelihood ratio'''
    L_H0 = log_likelihood(hist, 'H0', N_events, E_H0, FE2_H0, E_H1, FE2_H1)
    L_H1 = log_likelihood(hist, 'H1', N_events, E_H0, FE2_H0, E_H1, FE2_H1)
    LLR = L_H1 - L_H0
    return LLR


def plot_LLR_value_in_hist(N_events, bins, Eresolution, H, hist, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Function that fills the given histogram 'hist', with the LLR for 
    pseudo experiment data based on hypothesis 'H' '''
    h4 = pseudo_exp(N_events, bins, Eresolution, H, E_H0, FE2_H0, E_H1, FE2_H1)
    LLratio_Hdata = LLR(h4, N_events, E_H0, FE2_H0, E_H1, FE2_H1)
    hist.Fill(LLratio_Hdata)
    h4.Delete()
    return LLratio_Hdata


def pseudo_exp(N_events, bins, Eresolution, H, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Functin that creates pseudo experiments of N_events detections,
    based on the acceptence - rejection method'''
#    c3 = ROOT.TCanvas()
    h3 = ROOT.TH1D( 'h3','accept_reject, N=%s'%(N_events), bins, min(E_H0), max(E_H0))
    for i in range(N_events):
#        if i % (N_events/10) == 0:
#            print "%s/%s" %(i, N_events)
        E = accept_reject(E_H0, FE2_H0, E_H1, FE2_H1, H)
        E = E * ROOT.gRandom.Gaus(1., Eresolution/100) #Gauss(mean, sigma)          # <= wat doen met E<0?!
        h3.Fill(E)
#    h3.Draw()
    h3.SetXTitle("Energy (eV)")
    h3.SetYTitle("counts per bin")
#    c3.SaveAs("pseudo_event_Eresol%s_N%s_%s.png"%(Eresolution, N_events,H))
    return h3


def determine_p_value(h_h0, teststatistic_data):
    '''Funtion to determine the p-value of the data with given test statistic: 
    determine P(teststatistic > teststatistic_data | H0) by integration'''
    # integrate from tets_statistic to last bin
    axis      = h_h0.GetXaxis()
    bin_min   = axis.FindBin(teststatistic_data)
    bin_max   = h_h0.GetNbinsX() + 1
    integral  = h_h0.Integral(bin_min, bin_max)
    # subtract part of the first bin that was too much
    integral -= h_h0.GetBinContent(bin_min) * ( teststatistic_data - axis.GetBinLowEdge(bin_min) ) / axis.GetBinWidth(bin_min)
    # normalize
    p_value = integral/h_h0.GetEntries()
    return p_value


def determine_CL_value(h_h1, teststatistic_data):
    '''Funtion to determine the CL-value of the data with given test statistic: 
    determine P(teststatistic < teststatistic_data | H1) by integration'''
    # integrate from first bin to tets_statistic
    axis      = h_h1.GetXaxis()
    bin_min   = 1
    bin_max   = axis.FindBin(teststatistic_data)
    integral  = h_h1.Integral(bin_min, bin_max)
    # subtract part of the last bin that was too much
    integral -= h_h1.GetBinContent(bin_max) * ( (axis.GetBinLowEdge(bin_max) + h_h1.GetBinWidth(bin_max)) - teststatistic_data ) / axis.GetBinWidth(bin_min)
    # normalize
    CL_value = integral/h_h1.GetEntries()
    return CL_value


def plot_Hypothesis_test(h_h0, h_h1, savingname):
    cc = ROOT.TCanvas()

    #   Make up the plot
    h_h0.Draw()
    h_h0.SetLineColor(1)
    h_h1.Draw('SAME')
    h_h1.SetLineColor(3)
    #cc.BuildLegend()
    h_h1.SetTitle('Hypothesis Testing')
    h_h0.SetTitle('Hypothesis Testing')
    h_h0.SetXTitle("Test statistic log(lambda)")
    h_h0.SetYTitle("Probability density (counts)")

    cc.SaveAs(savingname)


def plot_line(canvas, hist, teststatistic, savingname):
    '''Function that plots the line of the test statistic of the pseudo experiment'''
    line = ROOT.TLine(teststatistic, 0, teststatistic, 3/5.*hist.GetMaximum())
    line.SetLineColor(2)                                                            # for shaded area: https://root.cern.ch/root/roottalk/roottalk98/0846.html
    line.SetLineWidth(3)
    line.Draw('SAME')
    canvas.Update()
    canvas.SaveAs(savingname)


def median(lst):
    '''Function to determine the median value of a given list 'lst' '''
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


# =================
# main & usage
# =================

def usage():
    print "Usage:  python  %s  <H0data-file>  <H1data-file> \n" % os.path.basename(sys.argv[0])

def main():
    
    ROOT.gStyle.SetOptStat(0)                                                   # do not show statistics box

    args = sys.argv[1:]
    if "-h" in args or "--help" in args or len(args) < 2:
        usage()
        sys.exit(2)

    # -- import flux data
    E_H0, FE2_H0 = np.loadtxt(sys.argv[1], delimiter = ',', unpack=True)
    E_H1, FE2_H1 = np.loadtxt(sys.argv[2], delimiter = ',', unpack=True)
#    draw_H0_H1(E_H0, FE2_H0, E_H1, FE2_H1)

    # -- normalize flux data (and plot)
    Binsize = E_H0[1]-E_H0[0]
    FE2_H1_norm = FE2_H1 / (sum(FE2_H1) * Binsize)
    FE2_H0_norm = FE2_H0 / (sum(FE2_H0) * Binsize)                              # Also divide by binsize to make sure that the interpolated function is normalized to 1
    #draw_H0_H1(E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)

    # -- check that the continuous function made with spectrum_func given FE2_H_norm data equals 1
#    f0_to_integrate = lambda e: spectrum_func(e, 'H0', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)
#    f1_to_integrate = lambda e: spectrum_func(e, 'H1', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)
#    integrant0, err1 = integrate.quad(f0_to_integrate, min(E_H0), max(E_H0))
#    integrant1, err1 = integrate.quad(f1_to_integrate, min(E_H0), max(E_H0))
#    print 'integrant0, integrant1', integrant0, integrant1

    # -- Create Pseudo Measurement for N_events according to H0, and H1, (plot, and save figure)
#    pseudo_exp(1000, 100, 0, 'H0', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)
#    pseudo_exp(1000, 100, 0, 'H1', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)


    # =================
    # -- Perform hypothesis testing:
    #    Create pseudo events based on H0 and H1. Determine, plot and store LLR's
    # =================

    #   Create needed canvas, empty histograms/lists, and set parameters
    ROOT.gStyle.SetTitleOffset(1.3, "y")

    N_events = []
    p_value = []
    CL_value = []

    for N in range(150):
        N = N+1
        if N % 5 == 0:
            print N


            N_events.append(N)      # for each pseudo experiment
            bins = 100              # amount of bins in histogram of pseudo experiment
            Eresolution = 0         # 0% error in Energy --> perfectly defined
    
            if 0 < N_events[-1] <= 10:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -5, 5)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -5, 5)      
            elif 10 < N_events[-1] <= 100:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -15, 15)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -15, 15)
            elif 100 < N_events[-1] <= 500:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -30, 30)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -30, 30)
            else:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -50, 50)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -50, 50)
            LLR_H0data = []
            LLR_H1data = []
    
            #   Check LLR for 'I_repetitions' pseudo events per hypothesis
            I_repetitions = 1000
            for i in range(I_repetitions):
                if i % 10 == 0:
                    print "i = %s/%s" %(i,I_repetitions)
                LLR_H0data.append(plot_LLR_value_in_hist(N_events[-1], bins, Eresolution, 'H0', h_h0, E_H0, FE2_H0_norm, E_H1, FE2_H1_norm))
                LLR_H1data.append(plot_LLR_value_in_hist(N_events[-1], bins, Eresolution, 'H1', h_h1, E_H0, FE2_H0_norm, E_H1, FE2_H1_norm))
    
            savingname = "H0_H1_data_Nevt_%i_Irep_%i.txt"%(N_events[-1], I_repetitions)
            with open(savingname, 'w') as f:
                for i in range(len(LLR_H0data)):
                    line = "%s,%s\n" %(LLR_H0data[i], LLR_H1data[i])
                    f.write(line)
    
            # plot and save Hypothesis testing
#            savingname = "Hypothesis_testing_Nevt_%i_Irep_%i.png"%(N_events[-1], I_repetitions)
#            plot_Hypothesis_test(h_h0, h_h1, savingname)

            # determine the median values of f(LLR|H1), and f(LLR|H0):
            median_H1 = median(LLR_H1data)
            median_H0 = median(LLR_H0data)

            # =================
            # -- P-value analysis
            # =================
            # determine expected p-value if H1 is true, print and plot
            p_value.append(determine_p_value(h_h0, median_H1))
            print 'expected p-value if H1 is true = ', p_value[-1]
            #plot_line(cc, h_h0, median_H1, "Hypothesis_testing_p_value_Nevt_%i_Irep_%i.png"%(N_events[-1], I_repetitions))

            # =================
            # -- CL(s+b)-value analysis
            # =================
            # determine expected p-value if H1 is true, print and plot
            CL_value.append(determine_CL_value(h_h0, median_H1))
            print 'expected CL-value if H1 can be excluded = ', CL_value[-1]
            #plot_line(cc, h_h1, median_H0, "Hypothesis_testing_CL_value_Nevt_%i_Irep_%i.png"%(N_events[-1], I_repetitions))
    
            h_h0.Delete()
            h_h1.Delete()

    print "Now write data to txt file"
    with open("N_p_CL_test.txt", 'w') as f:
        for i in range(len(N_events)):
            line = "%s,%s,%s\n" %(N_events[i], p_value[i], CL_value[i])
            f.write(line)

#    plt.plot(N_events, p_value)
#    plt.show()

#    plt.plot(N_events, CL_value)
#    plt.show()


if __name__ == "__main__":
    sys.exit(main())


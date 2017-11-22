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


def plot_LLR_value_in_hist(N_events, bins, H, hist, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Function that fills the given histogram 'hist', with the LLR for 
    pseudo experiment data based on hypothesis 'H' '''
    h4 = pseudo_exp(N_events, bins, H, E_H0, FE2_H0, E_H1, FE2_H1)
    LLratio_Hdata = LLR(h4, N_events, E_H0, FE2_H0, E_H1, FE2_H1)
    hist.Fill(LLratio_Hdata)
    h4.Delete()
    return LLratio_Hdata


def pseudo_exp(N_events, bins, H, E_H0, FE2_H0, E_H1, FE2_H1):
    '''Functin that creates pseudo experiments of N_events detections,
    based on the acceptence - rejection method'''
#    c3 = ROOT.TCanvas()
    h3 = ROOT.TH1D( 'h3','accept_reject, N=%s'%(N_events), bins, min(E_H0), max(E_H0))
    for i in range(N_events):
#        if i % (N_events/10) == 0:
#            print "%s/%s" %(i, N_events)
        E = accept_reject(E_H0, FE2_H0, E_H1, FE2_H1, H)
        h3.Fill(E)
    h3.Draw()
    h3.SetXTitle("Energy (eV)")
    h3.SetYTitle("counts per bin")
#    c3.SaveAs("check_pseudo_event_%s.png"%(H))
#    c3.Delete()
    return h3


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
    draw_H0_H1(E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)

    # -- check that the continuous function made with spectrum_func given FE2_H_norm data equals 1
#    f0_to_integrate = lambda e: spectrum_func(e, 'H0', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)
#    f1_to_integrate = lambda e: spectrum_func(e, 'H1', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)
#    integrant0, err1 = integrate.quad(f0_to_integrate, min(E_H0), max(E_H0))
#    integrant1, err1 = integrate.quad(f1_to_integrate, min(E_H0), max(E_H0))
#    print 'integrant0, integrant1', integrant0, integrant1

    # -- Create Pseudo Measurement for N_events according to H0, and H1, (plot, and save figure)
#    pseudo_exp(1000, 100, 'H0', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)
#    pseudo_exp(1000, 100, 'H1', E_H0, FE2_H0_norm, E_H1, FE2_H1_norm)


    # =================
    # -- Perform hypothesis testing:
    #    Create pseudo events based on H0 and H1. Determine, plot and store LLR's
    # =================

    #   Create needed canvas, empty histograms and lists
    cc = ROOT.TCanvas()
    ROOT.gStyle.SetTitleOffset(1.3, "y")

    N_events = 50         # for each pseudo experiment
    bins = 100              # amount of bins in histogram of pseudo experiment
   
    if 0 < N_events <= 10:
        h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -5, 5)
        h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -5, 5)      
    elif 10 < N_events <= 100:
        h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -15, 15)
        h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -15, 15)
    elif 100 < N_events <= 500:
        h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -30, 30)
        h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -30, 30)
    else:
        h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 100, -50, 50)
        h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 100, -50, 50)
    LLR_H0data = []
    LLR_H1data = []

    #   Check LLR for 'I_repetitions' pseudo events per hypothesis
    I_repetitions = 100
    for i in range(I_repetitions):
        print "i = %s/%s" %(i,I_repetitions)
        LLR_H0data.append(plot_LLR_value_in_hist(N_events, bins, 'H0', h_h0, E_H0, FE2_H0_norm, E_H1, FE2_H1_norm))
        LLR_H1data.append(plot_LLR_value_in_hist(N_events, bins, 'H1', h_h1, E_H0, FE2_H0_norm, E_H1, FE2_H1_norm))

    #   Make up the plot
    h_h0.Draw()
    h_h0.SetLineColor(1)
    h_h1.Draw('SAME')
    h_h1.SetLineColor(3)
    cc.BuildLegend()
    h_h1.SetTitle('Hypothesis Testing')
    h_h0.SetTitle('Hypothesis Testing')
    h_h0.SetXTitle("Test statistic log(lambda)")
    h_h0.SetYTitle("Probability density (counts)")

    savingname = "Hypothesis_testing_Nevt_%i_Irep_%i.png"%(N_events, I_repetitions)
    cc.SaveAs(savingname)

if __name__ == "__main__":
    sys.exit(main())


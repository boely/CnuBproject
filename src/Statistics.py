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
ROOT.gROOT.SetBatch(True)   # so it doesn't try to actively pop-up graphics on the screen <= needed for stoomboot!
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from bisect import bisect_left
from datetime import datetime

filename = os.path.splitext(os.path.basename(sys.argv[1]))[0]
basename = filename.split("H0H1")[0]

Elb = 5e20 # range of interest for the dip <- should implement in config file...
Eub = 5e22 # range of interest for the dip <- should implement in config file...

# =================
# Flux data
# =================

def spectrum_func(e, H, E_H, Flux_H0, Flux_H1 ):    # mean expected events
    '''Function that returns a function for the null hypothesis H = H0 
    (constant spectrum), and H = H1 hypothesis (dip in spectrum) '''

    # interpolate the separate datapoints
    f_H0 = interp1d(E_H, Flux_H0)
    f_H1 = interp1d(E_H, Flux_H1)

    if H == 'H0':
        return f_H0(e)
    elif H == 'H1':
        return f_H1(e)


def fill_hist( func, H, E_H, Flux_H0, Flux_H1):
    '''Function to plot a histogram according to a function func, given the 
    hypothesis H'''
    bins = 100
    h1 = ROOT.TH1D( 'h1','hypotheses', bins, min(E_H), max(E_H))

    Binwidth = (max(E_H)-min(E_H))/bins

    for be in range( 1, h1.GetNbinsX() + 1 ):
        e = h1.GetXaxis().GetBinCenter( be )
        b = h1.GetBin( be ) # global bin number
        z = func(e, H, E_H, Flux_H0, Flux_H1 )
        h1.SetBinContent( b, z )

    h1.SetXTitle("Energy (eV)")
    h1.SetYTitle("Flux (Arbitrary unit)")

    return h1


def draw_H0_H1(what, E_H, Flux_H0, Flux_H1):
    '''Function that draws H0 and H1 from imported data in a histogram'''
    c1 = ROOT.TCanvas()
    h2 = fill_hist( spectrum_func , 'H1', E_H, Flux_H0, Flux_H1)
    h2.Draw('LSAME')
    h2.SetTitle("Flux function")
    h1 = fill_hist( spectrum_func , 'H0', E_H, Flux_H0, Flux_H1)
    h1.SetLineStyle(3)
    h1.Draw('LSAME')
    h1.SetTitle("Flux function")
    c1.Update()

    h1.Write("Original_Flux_H0"+what)
    h2.Write("Original_Flux_H1"+what)
    c1.Write("Original_Flux_canvas"+what)

    h1.Delete()
    h2.Delete()

# =================
# Inverse Transform
# =================

def binary_search(a, x):
    '''Function that returns the position in list a closest to value x
    if x is exactly halfway two values, it takes the upper index,
    if x is below the lowest value, index 0 is returned
    if x is above the upper value, the last index is returned'''
    if x < a[0]:
        return 0
    elif x > a[-1]:
        return -1
    else:
        pos = bisect_left(a, x)  # find insertion position coming from left
        if float(a[pos]) == float(x):
            return pos
        elif a[pos] - x > x - a[pos-1]:
            return pos - 1
        else:
            return pos


def integr_spectrum_func(E_H, Flux_H):
    '''Function thet returns stepwise integrated Flux'''
    Int_Flux_H = [Flux_H[0]]
    for i in range(len(E_H)-1):
        i += 1
        Int_Flux_H.append(Int_Flux_H[-1]+Flux_H[i])
    return Int_Flux_H


def inverse_transform(E_H, Int_Flux_H0, Int_Flux_H1, H):
    '''Funtion that returns the energy corrsponding to the randomly picked
    flux value from the integrated flux
    Function can be improven by doing it not binwise but continuous'''
    if H == 'H1': 
        Int_Flux = Int_Flux_H1
    elif H == 'H0':
        Int_Flux = Int_Flux_H0

    # bounds of the fluxvalues to create random number in:
    Int_Flux_max = max(Int_Flux)

    rand_f = ROOT.gRandom.Rndm() * Int_Flux_max                                   # uniformly distributed number between 0 and Flux_max
    rand_e = E_H[binary_search(Int_Flux, rand_f)]

    return rand_e

# =================
# Likelihood
# =================

def log_likelihood(hist, H, N_events, E_H, Flux_H0, Flux_H1):
    '''Function that returns the log likelihood for a given hypothesis 'H' with
    respect to the given histogram 'hist' with your data. 
    NB: make sure Flux_H0 and Flux_H1, interpolated are normalized to 1'''
    loglik = 0.

    # -- loop over bins
    for i_bin in range( 1, hist.GetNbinsX() + 1 ):
        e = hist.GetBinCenter( i_bin )                                          # energy (centre of bin) 
        mu_bin = spectrum_func(e, H, E_H, Flux_H0, Flux_H1) * hist.GetBinWidth(1) * N_events     # theoretical amount of counts in bin
        Nevt_bin = hist.GetBinContent( i_bin )                                  # the amount of counts in bin found in the data of histogram h

        if ROOT.TMath.Poisson( Nevt_bin, mu_bin ) > 0:                          # <= check if&when this will happen
            loglik += ROOT.TMath.Log( ROOT.TMath.Poisson( Nevt_bin, mu_bin ) )
        else:
            print "Negative probability!!"
	    print "Found amount of events in bin: ",i_bin," = ", Nevt_bin
	    print "Expected amount of events in bin: ",i_bin, " = ", mu_bin
	    print "Poisson would give probability of: ", ROOT.TMath.Poisson( Nevt_bin, mu_bin )
    return loglik


def LLR(hist, N_events, E_H, Flux_H0, Flux_H1):
    '''Function that caluclates the likelihood ratio'''
    L_H0 = log_likelihood(hist, 'H0', N_events, E_H, Flux_H0, Flux_H1)
    L_H1 = log_likelihood(hist, 'H1', N_events, E_H, Flux_H0, Flux_H1)
    LLR = L_H1 - L_H0
    return LLR


def plot_LLR_value_in_hist(N_events, bins, Eresolution, H, hist, E_H, Flux_H0, Flux_H1, Int_Flux_H0, Int_Flux_H1):
    '''Function that fills the given histogram 'hist', with the LLR for 
    pseudo experiment data based on hypothesis 'H' '''
    h4 = pseudo_exp(N_events, bins, Eresolution, H, E_H, Int_Flux_H0, Int_Flux_H1)
    LLratio_Hdata = LLR(h4, N_events, E_H, Flux_H0, Flux_H1)
    hist.Fill(LLratio_Hdata)
    h4.Delete()
    return LLratio_Hdata


def pseudo_exp(N_events, bins, Eresolution, H, E_H, Int_Flux_H0, Int_Flux_H1):
    '''Function that creates pseudo experiments of N_events detections,
    based on the acceptence - rejection method'''
    h3 = ROOT.TH1D( 'h3','Pseudo Experiment, N=%s'%(N_events), int(bins), min(E_H), max(E_H))
    for i in range(int(N_events)):
#        if i % (N_events/10) == 0:
#            print "%s/%s" %(i, N_events)
        E = inverse_transform(E_H, Int_Flux_H0, Int_Flux_H1, H)
        Enieuw = E * ROOT.gRandom.Gaus(1., Eresolution/100.) #Gauss(mean, sigma)
        h3.Fill(Enieuw)
    h3.SetXTitle("Energy (eV)")
    h3.SetYTitle("Counts per bin")

#    h3.Write('Pseudo_exp_InvTr, %s, N=%s, Eres=%s'%(H, N_events, Eresolution))
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

    h_h0.Write(savingname+'H0')
    h_h1.Write(savingname+'H1')

    cc.Write(savingname)


def plot_line(canvas, hist, teststatistic, savingname):
    '''Function that plots the line of the test statistic of the pseudo experiment'''
    line = ROOT.TLine(teststatistic, 0, teststatistic, 3/5.*hist.GetMaximum())
    line.SetLineColor(2)                                                            # for shaded area: https://root.cern.ch/root/roottalk/roottalk98/0846.html
    line.SetLineWidth(3)
    line.Draw('SAME')
    canvas.Update()
    canvas.Write(savingname)


def median(lst):
    '''Function to determine the median value of a given list 'lst' '''
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


def Flux_for_Resolution(N_events, bins, ERES, E_H, Int_Flux_H0, Int_Flux_H1):
    '''Function that applies the energy resolution for H0 and H1 theory in a
    wide range (range of the <flux data file>), and returns the smeared
    flux in the energy range of interest: Elb < E < Eub'''
    Bins = int((max(E_H)-min(E_H)) * int(bins) / (Eub - Elb)) # scaling
    print "ERES = ", ERES
    hist0 = pseudo_exp(int(N_events), Bins, ERES, 'H0', E_H, Int_Flux_H0, Int_Flux_H1)
    hist0.Write('Pseudo_exp_InvTr, H0, N=%s, Eres=%s'%(N_events, ERES))
    hist1 = pseudo_exp(int(N_events), Bins, ERES, 'H1', E_H, Int_Flux_H0, Int_Flux_H1)
    hist1.Write('Pseudo_exp_InvTr, H1, N=%s, Eres=%s'%(N_events, ERES))

    E_H = []
    Flux_H0 = []
    Flux_H1 = []

    binwidth = hist0.GetBinWidth(1)

    for i_bin in range(hist0.GetNbinsX()):
        if (hist0.GetBinCenter(i_bin) + 0.5 * binwidth) >= Elb and (hist0.GetBinCenter(i_bin) - 0.5 * binwidth ) <= Eub:
            E_H.append(hist0.GetBinCenter(i_bin))
            Flux_H0.append(hist0.GetBinContent(i_bin))
    
#    i = 0
    for i_bin in range(hist1.GetNbinsX()):
        if (hist1.GetBinCenter(i_bin) + 0.5 * binwidth) >= Elb and (hist1.GetBinCenter(i_bin) - 0.5 * binwidth ) <= Eub:
#           if E_H[i] != hist1.GetBinCenter(i_bin):
#                print "EH[i] of H0 = ", E_H[i], " whereas E_H[i] of H1 = ", hist1.GetBinCenter(i_bin), "\n this is in bin i = ", i
            Flux_H1.append(hist1.GetBinContent(i_bin))
#            i += 1
    return np.asarray(E_H), np.asarray(Flux_H0), np.asarray(Flux_H1) # needed to be numpy arrays

# =================
# main & usage
# =================

def usage():
    print "Usage:  python  %s  <Fluxdata-file> \n" % os.path.basename(sys.argv[0])

def main():

    startTime = datetime.now()
    print basename

    f = ROOT.TFile(basename+"histos_RUNxx.root", "recreate")                     # create root file to save everything in   	  ################################# <========= RUN

    ROOT.gStyle.SetOptStat(0)                                                   # do not show statistics box
    ROOT.gStyle.SetTitleOffset(1.3, "y")                                        # spacing y-label

    args = sys.argv[1:]
    if "-h" in args or "--help" in args or len(args) < 1:
        usage()
        sys.exit(2)

    # =================
    # Flux data
    # =================

    # -- Import flux data
    E_H, Flux_H0, Flux_H1 = np.loadtxt(sys.argv[1], delimiter = '\t', unpack=True)
    draw_H0_H1("Flux", E_H, Flux_H0, Flux_H1)

    # -- Flux does not need to be normalized here yet (!).
    #    Inverse transform method only needs shape

    # -- Integrate flux data for inverse transform method
    Int_Flux_H0 = integr_spectrum_func(E_H, Flux_H0)
    Int_Flux_H1 = integr_spectrum_func(E_H, Flux_H1)
    draw_H0_H1("Int", E_H, Int_Flux_H0, Int_Flux_H1)

    # =================
    # Energy Resolution
    # =================
    bins = 100              # amount of bins in histogram of pseudo experiment (between Elb and Eub)

    N = int(1e6)                                                             				  ################################# <=========
    # -- Create & save ONE Smeared Pseudo Measurement with maaany events for new H0 and H1 shape including resolution
    if "ERES0_" in basename:
        print 'ERES0'
        E_H, Flux_H0, Flux_H1 = Flux_for_Resolution(N, bins, 0., E_H, Int_Flux_H0, Int_Flux_H1)
    elif "ERES30_" in basename:
        print 'ERES30'
        E_H, Flux_H0, Flux_H1 = Flux_for_Resolution(N, bins, 30., E_H, Int_Flux_H0, Int_Flux_H1)
    elif "ERES80_" in basename:
        print 'ERES80'
        E_H, Flux_H0, Flux_H1 = Flux_for_Resolution(N, bins, 80., E_H, Int_Flux_H0, Int_Flux_H1)
    else:
        print "Couldn't find que 'ERES0', 'ERES30', or 'ERES80' in basneme... CONTINUED WITH ERES0"
        E_H, Flux_H0, Flux_H1 = Flux_for_Resolution(N, bins, 0., E_H, Int_Flux_H0, Int_Flux_H1)

    # =================
    # Flux data with Eresolution applied
    # =================
    draw_H0_H1("FluxERES", E_H, Flux_H0, Flux_H1)

    # -- Redefine flux datapoints so that I have 10000 datapoints between
    #    min(E_H) == 5e20 (Elb), and max(E_H) == 5e22 (Eub) again.
    E = np.linspace(Elb, Eub, 10000) # <= from Neutrinoflux.py
    Flux_tmp_H0 = []
    Flux_tmp_H1 = []
    for e in E:
        Flux_tmp_H0.append(spectrum_func(e, 'H0', E_H, Flux_H0, Flux_H1))
        Flux_tmp_H1.append(spectrum_func(e, 'H1', E_H, Flux_H0, Flux_H1))
    E_H = E

    # -- Normalize flux data (and plot)
    #  - This is needed to easily determine the theoretically expected counts 
    #    in an Energy bin when comparing each pseudo event with the expected 
    #    flux data (H0&H1).
    f0_to_integrate = lambda e: spectrum_func(e, 'H0', E_H, Flux_tmp_H0, Flux_tmp_H1)
    f1_to_integrate = lambda e: spectrum_func(e, 'H1', E_H, Flux_tmp_H0, Flux_tmp_H1)
    integrant0, err1 = integrate.quad(f0_to_integrate, Elb, Eub)
    integrant1, err1 = integrate.quad(f1_to_integrate, Elb, Eub)

    Flux_H0_norm = np.asarray(Flux_tmp_H0) / integrant0
    Flux_H1_norm = np.asarray(Flux_tmp_H1) / integrant1
    draw_H0_H1("NormERES", E_H, Flux_H0_norm, Flux_H1_norm)

#    # -- Check that the continuous function made with spectrum_func given Flux_H_norm data equals 1
#    f0_to_integrate = lambda e: spectrum_func(e, 'H0', E_H, Flux_H0_norm, Flux_H1_norm)
#    f1_to_integrate = lambda e: spectrum_func(e, 'H1', E_H, Flux_H0_norm, Flux_H1_norm)
#    integrant0, err1 = integrate.quad(f0_to_integrate, Elb, Eub)
#    integrant1, err1 = integrate.quad(f1_to_integrate, Elb, Eub)
#    print 'integrant0, integrant1', integrant0, integrant1

    # -- Integrate flux data for inverse transform method
    Int_Flux_H0 = integr_spectrum_func(E_H, Flux_H0_norm)
    Int_Flux_H1 = integr_spectrum_func(E_H, Flux_H1_norm)
    draw_H0_H1("IntERES", E_H, Int_Flux_H0, Int_Flux_H1)

    # =================
    # -- Perform hypothesis testing:
    #    Create pseudo events based on H0 and H1. Determine, plot and store LLR's
    # =================

    #   Create needed canvas, empty histograms/lists, and set parameters
    ROOT.gStyle.SetTitleOffset(1.3, "y")

    N_events = []
    p_value = []
    CL_value = []

    for N in range(350):                                                                      ################################# <=========
        N = N+301
        if N % 50 == 0:              # option to analyze N stepwise: %1 analyses every N       ################################# <=========
            print "For %s events, create pseudo events." %N

            N_events.append(N)      # for each pseudo experiment
            Eresolution = 0         # Since it is already implemented in shape of H0 and H1.

            if 0 < N_events[-1] <= 10:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 1000, -5, 5)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 1000, -5, 5)      
            elif 10 < N_events[-1] <= 100:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 3000, -15, 15)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 3000, -15, 15)
            elif 100 < N_events[-1] <= 500:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 6000, -30, 30)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 6000, -30, 30)
            else:
                h_h0 = ROOT.TH1D( 'h_h0','LLR for H0', 10000, -50, 50)
                h_h1 = ROOT.TH1D( 'h_h1','LLR for H1', 10000, -50, 50)
            LLR_H0data = []
            LLR_H1data = []

            #   Check LLR for 'I_repetitions' pseudo events per hypothesis
            I_repetitions = int(1e4)                                                              ################################# <=========
            for i in range(I_repetitions):
                if i % 1000 == 0:
                    print "i = %s/%s" %(i,I_repetitions)
                LLR_H0data.append(plot_LLR_value_in_hist(N_events[-1], bins, Eresolution, 'H0', h_h0, E_H, Flux_H0_norm, Flux_H1_norm, Int_Flux_H0, Int_Flux_H1))
                LLR_H1data.append(plot_LLR_value_in_hist(N_events[-1], bins, Eresolution, 'H1', h_h1, E_H, Flux_H0_norm, Flux_H1_norm, Int_Flux_H0, Int_Flux_H1))

            savingname = "%s_H0_H1_data_Nevt_%i_Irep_%i.txt" %(basename, N_events[-1], I_repetitions)
#            with open(savingname, 'w') as txtfile:
#                for i in range(len(LLR_H0data)):
#                    line = "%s,%s\n" %(LLR_H0data[i], LLR_H1data[i])
#                    txtfile.write(line)

            # plot and save Hypothesis testing
            plot_Hypothesis_test(h_h0, h_h1, savingname)

            # determine the median values of f(LLR|H1), and f(LLR|H0):
            median_H1 = median(LLR_H1data)
            median_H0 = median(LLR_H0data)

            # =================
            # -- P-value analysis
            # =================
            # determine expected p-value if H1 is true, print and plot
            p_value.append(determine_p_value(h_h0, median_H1))
            print 'expected p-value if H1 is true = ', p_value[-1]
            #plot_line(cc, h_h0, median_H1, "Hypothesis_testing_p_value_Nevt_%i_Irep_%i"%(N_events[-1], I_repetitions))

            # =================
            # -- CL(s+b)-value analysis
            # =================
            # determine expected p-value if H1 is true, print and plot
            CL_value.append(determine_CL_value(h_h1, median_H0))
            print 'expected CL-value if H1 can be excluded = ', CL_value[-1]
            #plot_line(cc, h_h1, median_H0, "Hypothesis_testing_CL_value_Nevt_%i_Irep_%i"%(N_events[-1], I_repetitions))

            h_h0.Delete()
            h_h1.Delete()

    print "Now write data to txt file"
    with open("%s_RUNxx_N_p_CL_test.txt" %(basename), 'w') as txtfile:                                                     ################################# <========= RUN
        for i in range(len(N_events)):
            line = "%s,%s,%s\n" %(N_events[i], p_value[i], CL_value[i])
            txtfile.write(line)

    f.Close()

    print "DURATION OF STATISTICS.PY", datetime.now() - startTime

if __name__ == "__main__":
    sys.exit(main())


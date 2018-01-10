import sys
from matplotlib import pyplot as plt
import numpy as np
from array import array
import os
import ROOT
from math import *

filename = os.path.splitext(os.path.basename(sys.argv[1]))[0]
print filename
basename = filename.split("__")[0]
filename = sys.argv[1]
RUN = filename.split("\ERES")[0]
RUN = RUN.split("RUN")[1]

savename = "RUN"+RUN+"_"+basename
print savename

def make_plot(x_array, y_array, x_err, y_err, title, x_title, y_title, 
              draw, savingname, scale):
    '''Function to 'easily make a nice plot with errorbars'''

    c2 = ROOT.TCanvas()
    if scale == 'log':
        c2.SetLogy()
    c2.SetGrid()

    n = len(x_array)
    x = array('d', x_array)
    y = array('d', y_array)    
    x_er = array('d', x_err)
    y_er = array('d', y_err)

    gr = ROOT.TGraphErrors(n, x, y, x_er, y_er)

    gr.SetMarkerStyle(2)
    gr.SetMarkerColor(1)
    gr.SetTitle(title)
    gr.SetLineWidth(2)

    gr.Draw(draw)
    gr.GetXaxis().SetTitle(x_title)
    gr.GetXaxis().SetTitleOffset(1.5)
    gr.GetYaxis().SetTitle(y_title)
    gr.GetYaxis().SetTitleOffset(1.5)

    gr.Write(y_title+"_graph")
    c2.Write(y_title+"canvas")
    return gr, c2


def findN(y, p0, p1, p0_err, p1_err):
    N = (log(y)-p0)/float(p1)
    dN_dp0 = -1/p1
    dN_dp1 = (p0-log(y))/(p1**2)
    N_err = sqrt((dN_dp0 * p0_err)**2 + (dN_dp1 * p1_err)**2)
    return ceil(N), ceil(N_err)


def main():

    args = sys.argv[1:]
    if len(args) > 1:
        Irepetitions = float(sys.argv[2])
    else:
        Irepetitions = 10000.

    N, p_value, CL = np.loadtxt(sys.argv[1], delimiter = ',', unpack=True)
    CL_value = 1 - CL
    filename = os.path.splitext(os.path.basename(sys.argv[1]))[0]
    print filename
    
    # create root file to save everything in
    f = ROOT.TFile(savename+"_RESULTS.root", "recreate")

    N_err = [0.]*len(N)
    p_err = [0.]*len(N)
    for i in range(len(N)):
        p_err[i] = sqrt((p_value[i]*(1-p_value[i]))/Irepetitions)
    CL_err = [0.]*len(N)
    for i in range(len(N)):
        CL_err[i] = sqrt((CL_value[i]*(1-CL_value[i]))/Irepetitions)

    x_max = 1500.
    y_min = 4e-7

    ########################
    # P-value
    ########################
    p_gr, c2 = make_plot(N, p_value, N_err, p_err, "p-value", 
    "Number of neutrino detections", "p-value", "A*", "test_P.png", 'log')

    #RANGE
    axis = p_gr.GetXaxis()
    axis.SetLimits(0., x_max) # range on x-axis
    p_gr.GetHistogram().SetMinimum(y_min) # max on y-axis

    # FIT
    p_f = ROOT.TF1("p_f", "expo")#
    p_gr.Fit(p_f)
    p_gr.Draw("A*")
    p_fit = p_gr.GetFunction("p_f")
    p0 = p_fit.GetParameter(0)
    p0_err = p_fit.GetParError(0)
    p1 = p_fit.GetParameter(1)
    p1_err = p_fit.GetParError(1)

    N005, N005_err = findN(0.05, p0, p1, p0_err, p1_err)
    N3s, N3s_err = findN(2.7e-3, p0, p1, p0_err, p1_err)
    N5s, N5s_err = findN(5.7e-7, p0, p1, p0_err, p1_err)
    print "\nN for p=0.05 = %i +- %i" %(N005, N005_err)
    print "N for 3sigma = %i +- %i" %(N3s, N3s_err)
    print "N for 5sigma = %i +- %i\n" %(N5s, N5s_err)

    # p = 0.05
    line005 = ROOT.TLine(0, 0.05, x_max, 0.05)
    line005.SetLineColor(416)
    line005.SetLineWidth(2)
    line005.SetLineStyle(2)
    line005.Draw('SAME')
    t005 = ROOT.TText(x_max +20, 0.05, "p=0.05")
    t005.SetTextColor(417)
    t005.SetTextSize(0.035)
    t005.Draw('SAME')

    # 3 sigma
    line3s = ROOT.TLine(0, 2.7e-3, x_max, 2.7e-3)
    line3s.SetLineColor(416)
    line3s.SetLineWidth(2)
    line3s.SetLineStyle(2)
    line3s.Draw('SAME')
    t3s = ROOT.TLatex(x_max+20, 2.7e-3, "3#sigma")
    t3s.SetTextColor(417)
    t3s.SetTextSize(0.035)
    t3s.Draw('SAME')

    # 5 sigma
    line5s = ROOT.TLine(0, 5.7e-7, x_max, 5.7e-7)
    line5s.SetLineColor(416)
    line5s.SetLineWidth(2)
    line5s.SetLineStyle(2)
    line5s.Draw('SAME')
    t5s = ROOT.TLatex(x_max+20, 5.7e-7 ,"5#sigma")
    t5s.SetTextColor(417)
    t5s.SetTextSize(0.035)
    t5s.Draw('SAME') 
 
    c2.Update()
    p_gr.Write("p-value_graph_FIT")
    c2.Write("p-value_canvas_FIT")


    ########################
    # P-value
    ########################
    CL_gr, c2 = make_plot(N, CL_value, N_err, CL_err, "1-CL level", 
    "Number of neutrino detections", "1-CL value", "A*", "test_CL.png", '')

    #RANGE
    axis = CL_gr.GetXaxis()
    axis.SetLimits(0., x_max) # range on x-axis
    CL_gr.GetHistogram().SetMaximum(1) # max on y-axis

    # FIT
    CL_f = ROOT.TF1("CL_f", "1-exp([0]+[1]*x)")
    CL_f.SetParameters(p0, p1)
    CL_gr.Fit(CL_f)
    CL_gr.Draw("A*")


    CL_fit = CL_gr.GetFunction("CL_f")
    p0 = CL_fit.GetParameter(0)
    p0_err = p_fit.GetParError(0)
    p1 = CL_fit.GetParameter(1)
    p1_err = p_fit.GetParError(1)

    N095, N095_err = findN(1-0.95, p0, p1, p0_err, p1_err) 
    # 1-0.95 since y=1-exp(p0+p1*N) --> N = ((log(1-y)-p0)/p1)
    print "\nN for 1-CL=0.95 = %i +- %i" %(N095, N095_err)

    # cl = 0.95
    line095 = ROOT.TLine(0, 0.95, x_max, 0.95)
    line095.SetLineColor(416)
    line095.SetLineWidth(2)
    line095.SetLineStyle(2)
    line095.Draw('SAME')
    t095 = ROOT.TText(x_max -295, 0.95+0.01, "1-CL=0.95")
    t095.SetTextColor(417)
    t095.SetTextSize(0.035)
    t095.Draw('SAME')

    c2.Update()
    p_gr.Write("1-CL-value_graph_FIT")
    c2.Write("1-CL-value_canvas_FIT")

    f.Close()


if __name__ == "__main__":
    sys.exit(main())

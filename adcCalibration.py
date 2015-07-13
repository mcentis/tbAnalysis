from ROOT import *
import sys
from math import *

def landauFun(x, par):
    mpshift  = -0.22278298 # Landau maximum location
    mpc = par[1] - mpshift * par[0];
    return TMath.Landau(x[0],mpc,par[0]) / par[0]

if len(sys.argv) != 3:
    print '\tUsage: python adcCalibration runList directory'
    sys.exit(1)

runInfo = []

thickness = 100. # [um]
errThick = 1. # [um]
errAngle = radians(1.) # [rad]
ionizationMPV = 79. # e-h pairs / um [um-1]
ionizationMean = 98. # e-h pairs / um [um-1]

with open(sys.argv[1], 'r') as runInfoFile:
    for line in runInfoFile:
        runInfo.append(line.split())

#print runInfo

# ========================= stability of mean value of landau, it is not stable at all.

# land = TF1('land', landauFun, -50, 500000, 2)
# land.SetNpx(int(1e4))
# land.SetParameter(0, 3)
# land.SetParameter(1, 30)

# maxlist = [500, 1000, 5000, 10000, 50000, 100000]
# minp = -50

# for a in maxlist:
#     hist = TH1D('hist' + str(a), 'hist', int((a - minp) / 2.), minp, a)
#     for i in range((a - minp) * 5):
#         hist.Fill(land.GetRandom())

#     print 'max range %d\tentries %f\tmean %f' %(a, hist.GetEntries(), hist.GetMean())

calGrMPV = TGraphErrors()
calGrMPV.SetName('calGrMPV')
calGrMPV.SetTitle('Landau MPV')
calGrMPV.SetFillColor(kWhite)
calGrMPV.SetLineColor(kRed)
calGrMPV.SetMarkerColor(kRed)

calGrMean = TGraphErrors()
calGrMean.SetName('calGrMean')
calGrMean.SetName('Mean')
calGrMean.SetFillColor(kWhite)
calGrMean.SetLineColor(kBlue)
calGrMean.SetMarkerColor(kBlue)

for run in runInfo:
    print run
    fileName = sys.argv[2] + '/' + run[-1] + '.root'
    #print fileName
    inFile = TFile.Open(fileName)
    #print inFile
    hist = inFile.Get('signalDistrTimeCutDistCut')
    #print hist
    li = hist.GetListOfFunctions()
    #li.Print()
    func = li[0]#hist.GetFunction('gausLang')
    #print func
    print '\tMPV: %f +- %f' %(func.GetParameter(4), func.GetParError(4))
    angle = radians(float(run[2]))
    effThick = thickness / cos(angle)
    effThickErr = effThick * sqrt(pow(tan(angle) * errAngle, 2) + pow(errThick / thickness, 2))
    nPoint = calGrMPV.GetN()
    calGrMPV.SetPoint(nPoint, func.GetParameter(4), ionizationMPV * effThick)
    calGrMPV.SetPointError(nPoint, func.GetParError(4), ionizationMPV * effThickErr)

    hist = inFile.Get('signalDistrTimeCutDistCut_noisePeakSub')
    calGrMean.SetPoint(nPoint, hist.GetMean(), ionizationMean * effThick)
    calGrMean.SetPointError(nPoint, hist.GetMeanError(), ionizationMean * effThickErr)
    print '\tMean: %f +- %f\n' %(hist.GetMean(), hist.GetMeanError())

calFuncMPV = TF1('calFuncMPV', '[0] * x', 0, 100)

calCanMPV = TCanvas('calCanMPV', 'calCanMPV')
calGrMPV.Draw('ap')
#calGrMPV.Fit(calFuncMPV, 'r')
calGrMPV.GetXaxis().SetLimits(0, 100)
calGrMPV.GetYaxis().SetRangeUser(0, 15000)
calGrMPV.GetXaxis().SetTitle('Landau MPV [ADC]')
calGrMPV.GetYaxis().SetTitle('Expected ionization')

calFuncMean = TF1('calFuncMean', '[0] * x', 0, 200)

calCanMean = TCanvas('calCanMean', 'calCanMean')
calGrMean.Draw('ap')
#calGrMean.Fit(calFuncMean, 'r')
calGrMean.GetXaxis().SetLimits(0, 200)
calGrMean.GetYaxis().SetRangeUser(0, 25000)
calGrMean.GetXaxis().SetTitle('Landau Mean [ADC]')
calGrMean.GetYaxis().SetTitle('Expected ionization')

calGr = TMultiGraph()
calGr.Add(calGrMPV)
calGr.Add(calGrMean)

calFunc = TF1('calFuncMean', '[0] * x', 0, 200)
#calFunc = TF1('calFuncMean', 'pol1', 0, 200)

calCan = TCanvas('calCan', 'calCan')
calGr.Draw('ap')
calGr.Fit(calFunc, 'r')
calGr.GetXaxis().SetLimits(0, 110)
calGr.GetYaxis().SetRangeUser(0, 17000)
calGr.GetXaxis().SetTitle('Measured ionization [ADC]')
calGr.GetYaxis().SetTitle('Expected ionization')
calCan.BuildLegend()
calCan.Modified()
calCan.Update()

raw_input('...')

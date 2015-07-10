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

calGr = TGraphErrors()
calGr.SetName('calGr')

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
    print '\tMPV: %f +- %f\n' %(func.GetParameter(4), func.GetParError(4))
    angle = radians(float(run[2]))
    effThick = thickness / cos(angle)
    effThickErr = effThick * sqrt(pow(tan(angle) * errAngle, 2) + pow(errThick / thickness, 2))
    nPoint = calGr.GetN()
    calGr.SetPoint(nPoint, func.GetParameter(4), ionizationMPV * effThick)
    calGr.SetPointError(nPoint, func.GetParError(4), ionizationMPV * effThickErr)

calFunc = TF1('calFunc', '[0] * x', 0, 100)

calCan = TCanvas('calCan', 'calCan')
calGr.Draw('ap')
calGr.Fit(calFunc, 'r')
calGr.GetXaxis().SetLimits(0, 100)
calGr.GetYaxis().SetRangeUser(0, 15000)
calGr.GetXaxis().SetTitle('Landau MPV [ADC]')
calGr.GetYaxis().SetTitle('Expected ionization')

raw_input('...')

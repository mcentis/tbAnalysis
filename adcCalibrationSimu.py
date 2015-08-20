from ROOT import *
import sys
from math import *

if len(sys.argv) != 4:
    print '\tUsage: python adcCalibration runList directory confFile'
    sys.exit(1)

confDict = {}
with open(sys.argv[3], 'r') as confFile:
    for line in confFile:
        if len(line) == 0: continue
        sharpPos = line.find('#')
        line = line[:sharpPos]
        li = line.split()
        if len(li) == 0: continue
        if '=' in li: li.remove('=')
        #print li
        if len(li) != 2:
            print 'Error in reading conf file!!'
            print li
            continue
        confDict[li[0]] = li[1]

ADCtoe = float(confDict['ADCtoe'])
ADCtoeErr = float(confDict['ADCtoeErr'])

tCorr_p0 = float(confDict['tCorr_p0'])
tCorr_p1 = float(confDict['tCorr_p1'])
tCorr_p0Err = float(confDict['tCorr_p0Err'])
tCorr_p1Err = float(confDict['tCorr_p1Err'])
tCorr_p0p1Cov = float(confDict['tCorr_p0p1Cov'])

targetChipTemp = float(confDict['targetChipTemp'])
tempErr = float(confDict['tempErr'])

targetGain = tCorr_p0 + tCorr_p1 * targetChipTemp
targetGainErr = sqrt(pow(tCorr_p0Err, 2) + pow(tCorr_p1Err * targetChipTemp, 2) + 2 * tCorr_p0p1Cov * targetChipTemp)


runInfo = []

thickness = 100. # [um]
errThick = 1. # [um]
errAngle = radians(1.) # [rad]

ionMPVdict = {0:[7395, 47], 25:[8115, 49], 51:[12053, 80]}
ionMeanDict = {0:[9764, 66], 25:[10640, 65], 51:[15333, 84]}

with open(sys.argv[1], 'r') as runInfoFile:
    for line in runInfoFile:
        runInfo.append(line.split())

#print runInfo

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

ratioGr = TGraphErrors()
ratioGr.SetName('ratioGr')
ratioGr.SetMarkerStyle(20)

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
    angle = radians(float(run[2]))
    effThick = thickness / cos(angle)
    effThickErr = effThick * sqrt(pow(tan(angle) * errAngle, 2) + pow(errThick / thickness, 2))
    nPoint = calGrMPV.GetN()

    tempDistr = inFile.Get('tempDistr')
    tempTotErr = sqrt(pow(tempErr, 2) + pow(tempDistr.GetMeanError(), 2))
    gainMeas = tCorr_p0 + tCorr_p1 * tempDistr.GetMean()
    gainMeasErr = sqrt(pow(tCorr_p0Err, 2) + pow(tCorr_p1Err * tempDistr.GetMean(), 2) + pow(tCorr_p1 * tempTotErr, 2) + 2 * tCorr_p0p1Cov * tempDistr.GetMean())
    error = func.GetParameter(4) * sqrt(pow(targetGainErr / targetGain, 2) + pow(gainMeasErr / gainMeas, 2) + pow(func.GetParError(4) / func.GetParameter(4), 2))
    errMPV = error

    ion = ionMPVdict[int(run[2])]
    ionErr = ion[0] * sqrt(pow(ion[1] / ion[0], 2) + pow(effThickErr / effThick, 2))

    calGrMPV.SetPoint(nPoint, func.GetParameter(4), ion[0])
    calGrMPV.SetPointError(nPoint, error, ionErr)
    print '\tMPV: %f +- %f' %(func.GetParameter(4), error)

    hist = inFile.Get('signalDistrTimeCutDistCut_noisePeakSub')
    hist.GetXaxis().SetRangeUser(7, 511) # cut to improve mean determination

    error = hist.GetMean() * sqrt(pow(targetGainErr / targetGain, 2) + pow(gainMeasErr / gainMeas, 2) + pow(hist.GetMeanError() / hist.GetMean(), 2))

    ion = ionMeanDict[int(run[2])]
    ionErr = ion[0] * sqrt(pow(ion[1] / ion[0], 2) + pow(effThickErr / effThick, 2))

    calGrMean.SetPoint(nPoint, hist.GetMean(), ion[0])
    calGrMean.SetPointError(nPoint, error, ionErr)
    print '\tMean: %f +- %f\n' %(hist.GetMean(), error)
    
    ratio = func.GetParameter(4) / hist.GetMean()
    ratioErr = ratio * sqrt(pow(errMPV / func.GetParameter(4), 2) + pow(error / hist.GetMean(), 2))
    ratioGr.SetPoint(nPoint, angle, ratio)
    ratioGr.SetPointError(nPoint, errAngle, ratioErr)
    print '\tMPV / Mean: %f +- %f\n' %(ratio, ratioErr)

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

ratioCan = TCanvas('ratioCan', 'ratioCan')
ratioGr.Draw('ap')
ratioGr.GetXaxis().SetTitle('Angle [^{%circ}]')
ratioGr.GetYaxis().SetTitle('MPV / mean')

calGr = TMultiGraph()
calGr.Add(calGrMPV)
calGr.Add(calGrMean)

calFunc = TF1('calFuncMean', '[0] * x', 0, 200)
#calFunc = TF1('calFuncMean', 'pol1', 0, 200)

calCan = TCanvas('calCan', 'calCan')
calCan.SetGridx()
calCan.SetGridy()
calGr.Draw('ap')
calGr.Fit(calFunc, 'r')
calGr.GetXaxis().SetLimits(0, 110)
calGr.GetYaxis().SetRangeUser(0, 17000)
calGr.GetXaxis().SetTitle('Measured ionization [ADC]')
calGr.GetYaxis().SetTitle('Expected ionization [e^{-}]')
leg = calCan.BuildLegend()
leg.SetLineColor(kWhite)
calCan.Modified()
calCan.Update()

raw_input('...')

#!/usr/bin/env python
#
#  By Cameron Bravo <bravo@slac.stanford.edu>
#
#      Used to analyze histogramed data (HD) produced
#      by running makeHD on a .bin to produce a HD.root
#      from a icalScan run directory.
#
import numpy as np
import ROOT as r
from copy import deepcopy
from optparse import OptionParser
from DAQMap import layerToFeb

oPar = OptionParser()
oPar.add_option("-i", "--inDir", type="string", dest="inDir",
        default=".",help="Specify Input Filename", metavar="inDir")
oPar.add_option("-o", "--outfilename", type="string", dest="outfilename",
        default="anaIcalScan.root",help="Specify Output Filename", metavar="outfilename")
(options, args) = oPar.parse_args()

r.gROOT.SetBatch(True)

inDir = options.inDir

everything={}
pulseH = {}
pulseENC = {}
ph_g_dict={}
gain_g_dict={}
offset_g_dict={}

for layer in range(1,15):
    for module in range(4):
        index= layer+.1*module
        everything[index] = [{},{}] 
        for ical in range(0,105,5):
            pulseH[ical] = {}
            pulseENC[ical] = {}
            print 'Now Fitting Ical %i'%ical
            for cg in range(8):
                inFile = r.TFile( inDir )
                smData0_hh = deepcopy(getattr(inFile,"smData_%s_%s_hh" %(layer,module)))
                print layer, module
                for calCH in range(16):
                    for apv in range(1,5):
                        chan = 128*apv + 8*calCH + cg
                        if chan == 0: print ical, cg, calCH, apv
                        #Extract the pedestal from sample 0
                        scData0_h = deepcopy(smData0_hh.ProjectionY('scData0_ch%i_h'%(chan), chan+1, chan+1, "e"))
                        sampleMean = scData0_h.GetMean()
                        sampleNoise = scData0_h.GetRMS()
                        gaus_f = r.TF1('gaus_f','gaus', sampleMean-4*sampleNoise, sampleMean+4*sampleNoise)
                        gaus_f.SetParameter(0, 10.0)
                        gaus_f.SetParameter(1, sampleMean)
                        gaus_f.SetParameter(2, sampleNoise)
                        scData0_h.Fit(gaus_f, 'QR')
                        ped = gaus_f.GetParameter(1)
                        peak = gaus_f.GetParameter(1)
                        enc = gaus_f.GetParameter(2)
                        pulseH[ical][chan] = peak-ped
                        pulseENC[ical][chan] = enc
                        pass
                    pass
                inFile.Close()
                pass
            pass
            everything[index][0].update(pulseH)
            everything[index][1].update(pulseENC)
        
        ph_g = {}
        chList = []
        osList = []
        gainList = []

        for chan in range(128,640):
            phList = []
            encList = []
            icalList = []
            for ical in everything[index][0].keys():
                icalList.append(625.0 * ical)
                #icalList.append(float(ical))
                phList.append(everything[index][0][ical][chan])
                encList.append(everything[index][1][ical][chan])
                pass
            ph_g[chan] = r.TGraph( len(icalList), np.array(icalList), np.array(phList) )
            ph_g[chan].SetName('ph%i_g_%i'%(chan,index))
            ph_g[chan].SetTitle('Channel %i Pulse Heights;Ical [e^{-}];Pulse Height [ADC units]'%chan)
            ph_g[chan].SetMarkerStyle(3)
            p1_f = r.TF1('p1_f','[0]+[1]*x', 0.0, 5.0e4)
            fitOffset = p1_f.SetParameter(0,0.0)
            fitGain = p1_f.SetParameter(1,0.0)
            ph_g[chan].Fit(p1_f,'QR')
            fitOffset = p1_f.GetParameter(0)
            fitGain = p1_f.GetParameter(1)
            chList.append(float(chan))
            gainList.append(fitGain)
            osList.append(fitOffset)
            pass

        gain_g = r.TGraph( len(chList), np.array(chList), np.array(gainList) )
        gain_g.SetName('gain_g')
        gain_g.SetTitle('Gains by Channel;Channel;Gain [ADC units / e^{-}]')
        gain_g.SetMarkerStyle(3)
        gain_g_dict[index]=gain_g

        offset_g = r.TGraph( len(chList), np.array(chList), np.array(osList) )
        offset_g.SetName('offset_g')
        offset_g.SetTitle('Offsets by Channel;Channel;Offset [ADC units]')
        offset_g.SetMarkerStyle(3)
        offset_g_dict[index]=offset_g

        ph_g_dict[index]={}
        ph_g_dict[index].update(ph_g)
    pass
pass


outFile = r.TFile(options.outfilename,"RECREATE")
outFile.cd()

for gain in gain_g_dict:
    gain.Write()
    pass
for offset in offset_g_dict:
    offset.Write()
    pass
for ph in ph_g_dict:
    for chan in range(128,640): ph_g[chan].Write()
    pass
outFile.Close()

exit()

inFile = r.TFile(options.infilename)
smData_hh = []
scData_h = []

#Extract Data
print 'Extracting Data'
smDataSum_hh = deepcopy(inFile.smDataSum_hh)
scDataSum_h = []
for ss in range(6):
    smData_hh.append( deepcopy(getattr(inFile,'smData%i_hh'%ss)) )
    scData_h.append([])
    for cc in range(640):
        scData_h[ss].append(smData_hh[ss].ProjectionY('scData%i_ch%i_h'%(ss,cc), cc+1, cc+1, "e"))
        if ss==0: scDataSum_h.append(smDataSum_hh.ProjectionY( 'scDataSum%i_h'%cc, cc+1, cc+1, "e" ))
        pass
    pass
pulG = []
pulGe = []
pulB = []
pulBe = []
tbin = []
zeros = []
for ss in range(6): 
    tbin.append(float(ss))
    zeros.append(0.0)
    pass
#inFile.Close()

#Analyze Data
print 'Analyzing Data'
outFilename = (options.infilename[:-7] + 'anaHD.root')
print outFilename
outFile = r.TFile(outFilename,"RECREATE")
chList = []
pedsSum = []
fitPedsSum = []
peds = []
fitPeds = []
s2sPedDiff = []
noiseSum = []
fitNoiseSum = []
noiseRatioSum = []
noises = []
fitNoises = []
fitXXoNDF = []
for ss in range(6): 
    peds.append([])
    fitPeds.append([])
    s2sPedDiff.append([])
    noises.append([])
    fitNoises.append([])
    pass

for cc in range(640): 
    chList.append(float(cc))
    sumMean = scDataSum_h[cc].GetMean()
    sumNoise = scDataSum_h[cc].GetRMS()
    gaus_f = r.TF1('gaus_f','gaus', sumMean-4*sumNoise, sumMean+4*sumNoise)
    gaus_f.SetParameter(0, 10.0)
    gaus_f.SetParameter(1, sumMean)
    gaus_f.SetParameter(2, sumNoise)
    scDataSum_h[cc].Fit(gaus_f, 'QR')
    fitSumMean = gaus_f.GetParameter(1)
    fitSumNoise = gaus_f.GetParameter(2)
    fitXX = gaus_f.GetChisquare()
    fitNDF = gaus_f.GetNDF()
    for ss in range(6):
        sampleMean = scData_h[ss][cc].GetMean()
        sampleNoise = scData_h[ss][cc].GetRMS()
        gaus_f = r.TF1('gaus_f','gaus', sampleMean-4*sampleNoise, sampleMean+4*sampleNoise)
        gaus_f.SetParameter(0, 10.0)
        gaus_f.SetParameter(1, sampleMean)
        gaus_f.SetParameter(2, sampleNoise)
        scData_h[ss][cc].Fit(gaus_f, 'QR')
        fitMean = gaus_f.GetParameter(1)
        fitNoise = gaus_f.GetParameter(2)
        if cc == 540: 
            pulG.append(fitMean)
            pulGe.append(fitNoise)
            pass
        if cc == 570: 
            pulB.append(fitMean)
            pulBe.append(fitNoise)
            pass
        peds[ss].append(sampleMean)
        fitPeds[ss].append(fitMean)
        s2sPedDiff[ss].append(fitMean-fitSumMean)
        noises[ss].append(sampleNoise)
        fitNoises[ss].append(fitNoise)
        pass
    pedsSum.append(sumMean)
    noiseSum.append(sumNoise)
    #noiseRatioSum.append(sumNoise/fitSumNoise)
    fitPedsSum.append(fitSumMean)
    fitNoiseSum.append(fitSumNoise)
    #fitXXoNDF.append(fitXX/fitNDF)
    pass

pedsSum_g = r.TGraph(len(chList), np.array(chList), np.array(pedsSum))
pedsSum_g.SetName('pedsSum_g')
pedsSum_g.SetTitle('Pedestals;Physical Channel # ;Mean [ADC units]')
fitPedsSum_g = r.TGraph(len(chList), np.array(chList), np.array(fitPedsSum))
fitPedsSum_g.SetName('fitPedsSum_g')
fitPedsSum_g.SetTitle('Fit Pedestals;Physical Channel # ;Mean [ADC units]')

noiseSum_g = r.TGraph(len(chList), np.array(chList), np.array(noiseSum))
noiseSum_g.SetName('noiseSum_g')
noiseSum_g.SetTitle('Noise;Physical Channel # ;RMS [ADC units]')

fitNoiseSum_g = r.TGraph(len(chList), np.array(chList), np.array(fitNoiseSum))
fitNoiseSum_g.SetName('fitNoiseSum_g')
fitNoiseSum_g.SetTitle('Fit Noise;Physical Channel # ;RMS [ADC units]')

fitPvNSum_g = r.TGraph(len(fitPedsSum), np.array(fitPedsSum), np.array(fitNoiseSum))
fitPvNSum_g.SetName('fitPvNSum_g')
fitPvNSum_g.SetTitle('Fit Noise;Fit Pedestal [ADC units];Fit ENC [ADC units]')
#noiseRatioSum_g = r.TGraph(len(chList), np.array(chList), np.array(noiseRatioSum))
#noiseRatioSum_g.SetName('noiseRatioSum_g')
#noiseRatioSum_g.SetTitle('Fit Noise over Noise;Physical Channel # ;Noise Ratio')

#fitXXoNDF_g = r.TGraph(len(chList), np.array(chList), np.array(fitXXoNDF))
#fitXXoNDF_g.SetName('fitXXoNDF_g')
#fitXXoNDF_g.SetTitle('Fit #chi^{2} over NDF;Physical Channel # ;Goodness of Fit')

pulG_ge = r.TGraphErrors(len(tbin), np.array(tbin), np.array(pulG), np.array(zeros), np.array(pulGe))
pulG_ge.SetName('pulG_ge')
pulG_ge.SetTitle('Good Channel Pulse;Sample;Pulse Height [ADC units]')
pulB_ge = r.TGraphErrors(len(tbin), np.array(tbin), np.array(pulB), np.array(zeros), np.array(pulBe))
pulB_ge.SetName('pulB_ge')
pulB_ge.SetTitle('Bad Channel Pulse;Sample;Pulse Height [ADC units]')

peds_g = []
noises_g = []
fitPeds_g = []
s2sPedDiff_g = []
fitNoises_g = []
for ss in range(6):
    peds_g.append(r.TGraph(len(chList), np.array(chList), np.array(peds[ss])))
    peds_g[ss].SetName('peds%i_g'%ss)
    peds_g[ss].SetTitle('Sample %i Pedestal;Physical Channel # ;Mean [ADC units]'%ss)
    fitPeds_g.append(r.TGraph(len(chList), np.array(chList), np.array(fitPeds[ss])))
    fitPeds_g[ss].SetName('fitPeds%i_g'%ss)
    fitPeds_g[ss].SetTitle('Sample %i Fit Pedestals;Physical Channel # ;Fit Mean [ADC units]'%ss)
    s2sPedDiff_g.append(r.TGraph(len(chList), np.array(chList), np.array(s2sPedDiff[ss])))
    s2sPedDiff_g[ss].SetName('s2sPedDiff%i_g'%ss)
    s2sPedDiff_g[ss].SetTitle('Sample-to-Sample %i Pedestal Difference;Physical Channel # ;Sample Pedestal minus Sum Pedestal [ADC units]'%ss)
    noises_g.append(r.TGraph(len(chList), np.array(chList), np.array(noises[ss])))
    noises_g[ss].SetName('noises%i_g'%ss)
    noises_g[ss].SetTitle('Sample %i Noise;Physical Channel # ;Noise [ADC units]'%ss)
    fitNoises_g.append(r.TGraph(len(chList), np.array(chList), np.array(fitNoises[ss])))
    fitNoises_g[ss].SetName('fitNoises%i_g'%ss)
    fitNoises_g[ss].SetTitle('Sample %i Fit Noise;Physical Channel # ;Fit Noise [ADC units]'%ss)
    pass


outFile.cd()
smDataSum_hh.Write()
pedsSum_g.Write()
noiseSum_g.Write()
fitPedsSum_g.Write()
fitPvNSum_g.Write()
fitNoiseSum_g.Write()
#noiseRatioSum_g.Write()
#fitXXoNDF_g.Write()
for ss in range(6): 
    smData_hh[ss].Write()
    peds_g[ss].Write()
    noises_g[ss].Write()
    s2sPedDiff_g[ss].Write()
    pass
pulG_ge.Write()
pulB_ge.Write()

outFile.Close()
inFile.Close()







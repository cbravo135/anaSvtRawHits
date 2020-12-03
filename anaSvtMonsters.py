#!/bin/env python
import ROOT as r
from DAQMap import *

r.gSystem.Load("libevent")

myF = r.TFile("hps_010494/hpssvt_010494_00006_raw.root")

myT = myF.HPS_Event

nHits_h = {}
nHits_hh = {}
histName = "nHits_hw_hh"
nHits_hh[histName] = r.TH2D(histName,"Number of Hit Strips;FEB;Hybrid",10,-0.5,9.5,4,-0.5,3.5)
histName = "nHits_sw_hh"
nHits_hh[histName] = r.TH2D(histName,"Number of Hit Strips;Layer;Module",14,-0.5,13.5,4,-0.5,3.5)
histName = "nHits_h"
nHits_h[histName] = r.TH1D(histName,"Number of Hits in Event;Number of Hits;Events",12000,0.5,12000.5)
histName = "nMonFebs_h"
nHits_h[histName] = r.TH1D(histName,"Number of FEBs with Monster Hybrid;Number of FEBs;Events",10,0.5,10.5)
histName = "monHyb_hh"
nHits_hh[histName] = r.TH2D(histName,"Monster Hybrids;FEB;Hybrid",10,-0.5,9.5,4,-0.5,3.5)
histName = "monFebs_h"
nHits_h[histName] = r.TH1D(histName,"FEBs with Monster Hybrids;FEB;Events",10,-0.5,9.5)
histName = "lastMon_h"
nHits_h[histName] = r.TH1D(histName,"Time to Previous Event;Time [ns];Events",1000000,0.5,1000000.5)
histName = "lastMonFront_h"
nHits_h[histName] = r.TH1D(histName,"Time to Previous Event (FEB<4);Time [ns];Events",1000000,0.5,1000000.5)
histName = "lastMonBack_h"
nHits_h[histName] = r.TH1D(histName,"Time to Previous Event (FEB>3);Time [ns];Events",1000000,0.5,1000000.5)
histName = "nextMon_h"
nHits_h[histName] = r.TH1D(histName,"Time to Previous Trigger;Time [ns];Monster Events",1000000,0.5,1000000.5)
histName = "trigTime_h"
nHits_h[histName] = r.TH1D(histName,"Time to Previous Trigger;Time [ns];Monster Events",1000000,0.5,10000000.5)
for feb in range(10):
    histName = 'nHits%i_h'%(feb)
    nHits_h[histName] = r.TH1D(histName,"FEB %i;Number of Hits;Events"%feb,2560,0.5,2560.5)
    histName = 'nMonHyb%i_h'%(feb)
    nHits_h[histName] = r.TH1D(histName,"FEB %i;Number of Monster Hybrids;Events"%feb,4,0.5,4.5)
    for hyb in range(4):
        histName = 'nHits%i_%i_h'%(feb,hyb)
        nHits_h[histName] = r.TH1D(histName,"FEB %i Hybrid %i;Number of Hits;Events"%(feb,hyb),640,0.5,640.5)
        histName = 'monOcc%i_%i_h'%(feb,hyb)
        nHits_h[histName] = r.TH1D(histName,"FEB %i Hybrid %i;Physical Strip;Hits in Monster Events"%(feb,hyb),640,-0.5,639.5)
        pass
    pass

#for lay in range(1,15):
#    for mod in range(4):
#        if lay < 9 and mod > 1: continue
#        SWfloat = lay+0.1*mod
#        HWfloat = layerToFeb[SWfloat]
#        feb = round(HWfloat)
#        hyb = int((HWfloat - float(feb))*10.0)
#        histName = 'nHits%i_%i_h'%(feb,hyb)
#        print "SW: %f\tHW: %f\tLayer: %i\tModule: %i\tFEB: %i\tHybrid: %i"%(SWfloat,HWfloat,lay,mod,feb,hyb)
#        print int((HWfloat - float(feb))*10.0)
#        print int(round((HWfloat - float(feb))*10.0))
#        print (HWfloat - float(feb))*10.0
#        pass
#    pass
#exit(0)

EN = 0
evTime = 0
lastEvTime = 0
N = myT.GetEntries()
lastEvMon = False
print "Running on %i events"%N
for ev in myT:
    EN = EN + 1
    if EN%(N/100) == 0: print 'Event: %i'%EN
    rawHits = ev.SVTRawTrackerHits
    lastlastEvTime = lastEvTime
    lastEvTime = evTime
    evHead = ev.EventHeader
    evTime = evHead.getEventTime()
    febNhits = {}
    hybNhits = {}
    hybHits = {}
    for feb in range(10): 
        febNhits[feb] = 0
        hybNhits[feb] = {}
        hybHits[feb] = {}
        for hyb in range(4):
            hybNhits[feb][hyb] = 0
            hybHits[feb][hyb] = []
            pass
        pass
    nHits = 0
    for hit in rawHits:
        lay = hit.getLayer()
        mod = hit.getModule()
        strip = hit.getStrip()
        adcs = hit.getADCs()
        HWfloat = layerToFeb[lay+0.1*mod]
        feb = round(HWfloat)
        hyb = int(round((HWfloat - feb)*10))
        nHits += 1
        febNhits[feb] += 1
        hybNhits[feb][hyb] += 1
        hybHits[feb][hyb].append(strip)
        nHits_hh["nHits_hw_hh"].Fill(feb,hyb)
        nHits_hh["nHits_sw_hh"].Fill(lay,mod)
        pass
    nHits_h["nHits_h"].Fill(nHits)
    nMonFebs = 0
    for feb in range(10):
        histName = "nHits%i_h"%feb
        nHits_h[histName].Fill(febNhits[feb])
        nMonHyb = 0
        for hyb in range(4):
            histName = "nHits%i_%i_h"%(feb,hyb)
            nHits_h[histName].Fill(hybNhits[feb][hyb])
            if hybNhits[feb][hyb] > 100: 
                nMonHyb += 1
                nHits_hh["monHyb_hh"].Fill(feb,hyb)
                for ss in hybHits[feb][hyb]: nHits_h["monOcc%i_%i_h"%(feb,hyb)].Fill(ss)
                pass
            pass
        nHits_h['nMonHyb%i_h'%(feb)].Fill(nMonHyb)
        if nMonHyb > 0: 
            nMonFebs += 1
            nHits_h["monFebs_h"].Fill(feb)
            if feb < 4: nHits_h["lastMonFront_h"].Fill(evTime - lastEvTime)
            else: nHits_h["lastMonBack_h"].Fill(evTime - lastEvTime)
        pass
    nHits_h["trigTime_h"].Fill(evTime - lastEvTime)
    if lastEvMon: nHits_h["nextMon_h"].Fill(evTime - lastEvTime)
    if nMonFebs > 0: 
        lastEvMon = True
        nHits_h["lastMon_h"].Fill(evTime - lastEvTime)
        nHits_h["nMonFebs_h"].Fill(nMonFebs)
    else:
        lastEvMon = False
    pass

print adcs

outF = r.TFile("monsters/hps_010494_00006_anaSvtMonsters.root","RECREATE")
fDirs = {}
fDirs['overview'] = outF.mkdir('overview')
fDirs['febs'] = outF.mkdir('febs')
fDirs['hyb'] = outF.mkdir('hyb')
fDirs['hybOcc'] = outF.mkdir('hybOcc')
fDirs['overview'].cd()
nHits_hh["nHits_hw_hh"].Write()
nHits_hh["nHits_sw_hh"].Write()
nHits_h["nHits_h"].Write()
nHits_h["trigTime_h"].Write()
nHits_h["lastMon_h"].Write()
nHits_h["lastMonFront_h"].Write()
nHits_h["lastMonBack_h"].Write()
nHits_h["nextMon_h"].Write()
nHits_hh["monHyb_hh"].Write()
nHits_h["monFebs_h"].Write()
fDirs['febs'].cd()
nHits_h["nMonFebs_h"].Write()
for feb in range(10):
    fDirs['febs'].cd()
    histName = 'nHits%i_h'%(feb)
    nHits_h[histName].Write()
    nHits_h['nMonHyb%i_h'%(feb)].Write()
    for hyb in range(4):
        fDirs['hyb'].cd()
        histName = 'nHits%i_%i_h'%(feb,hyb)
        nHits_h[histName].Write()
        fDirs['hybOcc'].cd()
        histName = 'monOcc%i_%i_h'%(feb,hyb)
        nHits_h[histName].Write()
        pass
    pass
outF.Close()

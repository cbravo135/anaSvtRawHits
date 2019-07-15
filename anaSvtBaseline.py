#!/bin/env python
import ROOT as r

r.gSystem.Load("libevent")

myF = r.TFile("trees/hpssvt_009709_baseline.root")

myT = myF.HPS_Event


smData_hh = {}
locData_hh = r.TH2D('locData_hh','locData_hh',14,0.5,14.5,4,-0.5,3.5)
for lay in range(1,15):
    for mod in range(4):
        pName = 'smData_%i_%i_hh'%(lay,mod)
        smData_hh[pName] = r.TH2D(pName,pName,640,0,640,8000,0,16000)
        pass
    pass
               
EN = 0
N = myT.GetEntries()
for ev in myT:
    EN = EN + 1
    if EN%(N/100) == 0: print 'Event: %i'%EN
    rawHits = ev.SVTRawTrackerHits
    for hit in rawHits:
        lay = hit.getLayer() 
        mod = hit.getModule() 
        strip = hit.getStrip()
        adcs = hit.getADCs()
        pName = 'smData_%i_%i_hh'%(lay,mod)
        locData_hh.Fill(float(lay),float(mod))
        for sample in range(6): smData_hh[pName].Fill(float(strip),float(adcs[sample]))
        pass
    pass
outF = r.TFile("hh/hpssvt_009612_hh.root","RECREATE")
outF.cd()
locData_hh.Write()
for lay in range(1,15):
    for mod in range(4):
        pName = 'smData_%i_%i_hh'%(lay,mod)
        smData_hh[pName].Write()
        pass
    pass
outF.Close()

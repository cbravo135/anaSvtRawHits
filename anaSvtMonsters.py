#!/bin/env python
import ROOT as r
from DAQMap import *

r.gSystem.Load("libevent")

myF = r.TFile("run2019/hps_run9600_RawSvtHitsTree.root")

myT = myF.HPS_Event

for lay in range(1,15):
    for mod in range(4):
        if lay < 9 and mod > 1: continue
        HWfloat = layerToFeb[lay+0.1*mod]
        feb = round(HWfloat)
        hyb = (HWfloat - feb)*10
        print "Layer: %i\tModule: %i\tFEB: %i\tHybrid: %i"%(lay,mod,feb,hyb)
        pass
    pass


outF = r.TFile("monsters/hps_009600_SvtMonsters.root","RECREATE")
outF.cd()
outF.Close()

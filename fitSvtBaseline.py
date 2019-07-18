#!/usr/bin/env python
#
#  By Cameron Bravo <bravo@slac.stanford.edu>, modified by Aster Taylor <ataylor@slac.stanford.edu>
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
means={}
rmss={}
mean_g_dict={}
rms_g_dict={}

for layer in xrange(1,15):
    for module in xrange(4):
        if layer < 9 and module > 1: continue
        index= layer+.1*module
        everything[index] = [{},{}] 
        for chan in xrange(1,641):
            if layer<5 and chan>512: continue
            means[chan] = {}
            rmss[chan] = {}
            inFile = r.TFile( inDir )
            smData0_hh = deepcopy(getattr(inFile,"smData_%s_%s_hh" %(layer,module)))
            print layer, module
            scData0_h = deepcopy(smData0_hh.ProjectionY('scData0_ch%i_h'%(chan), chan+1, chan+1, "e"))
            sampleMean = scData0_h.GetMean()
            sampleNoise = scData0_h.GetRMS()
            gaus_f = r.TF1('gaus_f','gaus', sampleMean-4*sampleNoise, sampleMean+4*sampleNoise)
            gaus_f.SetParameter(0, 10.0)
            gaus_f.SetParameter(1, sampleMean)
            gaus_f.SetParameter(2, sampleNoise)
            scData0_h.Fit(gaus_f, 'QR')
            mean = gaus_f.GetParameter(1)
            rms = gaus_f.GetParameter(2)
            means[chan] = mean
            rmss[chan] = rms
            inFile.Close()
            pass
        everything[index][0].update(means)
        everything[index][1].update(rmss)

        chList=[]
        meanList=[]
        rmsList=[]
        null=[]
        for chan in xrange(1,641):
            if layer<5 and chan>512: continue
            chList.append(float(chan))
            meanList.append(float(everything[index][0][chan]))
            rmsList.append(float(everything[index][1][chan]))
            null.append(0.0)
            pass

        mean_g = r.TGraphErrors( len(chList), np.array(chList), np.array(meanList), np.array(null), np.array(rmsList) )
        mean_g.SetName('mean_%i_g'%index)
        mean_g.SetTitle('Means vs Channel for %i;Channel;Mean'%index)
        mean_g.SetMarkerStyle(3)
        mean_g_dict[index]=mean_g

        rms_g = r.TGraph( len(chList), np.array(chList), np.array(rmsList) )
        rms_g.SetName('rms_%i_g'%index)
        rms_g.SetTitle('rmss vs Channel for %i;Channel;rms'%index)
        rms_g.SetMarkerStyle(3)
        rms_g_dict[index]=rms_g
        pass
    pass


outFile = r.TFile(options.outfilename,"RECREATE")
outFile.cd()

for avg in mean_g_dict:
    mean_g_dict[avg].Write()
    pass
for sig in rms_g_dict:
    rms_g_dict[sig].Write()
    pass
outFile.Close()

exit

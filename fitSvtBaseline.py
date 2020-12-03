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
oPar.add_option("-i", "--inFile", type="string", dest="inFile",
        default="hh/hpssvt_009617_hh.root",help="Specify Input Filename", metavar="inFile")
oPar.add_option("-o", "--outfilename", type="string", dest="outfilename",
        default="fitSvtBaseline.root",help="Specify Output Filename", metavar="outfilename")
oPar.add_option("-n", "--inRun", type="int", dest="inRun", default=0,
        help="Run ID", metavar="inRun")
(options, args) = oPar.parse_args()

r.gROOT.SetBatch(True)

if options.inRun==0: 
    infilename = options.inFile
else: 
    infilename = "hh/hpssvt_%06d_hh.root"%options.inRun

everything={}
means={}
rmss={}
mean_g_dict={}
rms_g_dict={}
nIndices=[]
#for layer in xrange(1,15):
#    for module in xrange(4):
#        if layer < 9 and module > 1: continue
#        sw_f = layer+.1*module
#        l2f=layerToFeb[sw_f]
#        feb=round(layerToFeb[sw_f])
#        hybrid=int((layerToFeb[sw_f]-feb+0.005)*10.0)
#        nsw_f=feb+0.1*hybrid
#        hw_s='F%iH%i'%(feb,hybrid)
#        print 'sw_f: %f l2f: %f feb: %i hyb: %i nsw_f: %f hw_s: %s'%(sw_f, l2f, feb, hybrid, nsw_f, hw_s)
#        pass
#    pass
#
#exit(0)

inFile = r.TFile( infilename )
smData_hh = {}
for layer in xrange(1,15):
    for module in xrange(4):
        if layer < 9 and module > 1: continue
        sw_f = layer+.1*module
        feb=round(layerToFeb[sw_f])
        hybrid=int((layerToFeb[sw_f]-feb+0.005)*10.0)
        nsw_f=feb+0.1*hybrid
        hw_s='F%iH%i'%(feb,hybrid)
        nIndices.append(hw_s)
        if hasattr(inFile, "smData_%s_%s_hh" %(layer,module))==False: continue
        everything[sw_f] = [{},{}] 
        mean = {}
        rms = {}
        smData_hh[hw_s] = deepcopy(getattr(inFile,"smData_%s_%s_hh" %(layer,module)))
        smData_hh[hw_s].SetName("smData_%s_hh" %(hw_s))
        for chan in xrange(640):
            print "Channel:", chan
            if layer<5 and chan>512: continue
            means[chan] = {}
            rmss[chan] = {}
            print "Opening File"
            print layer, module, sw_f
            scData_h = deepcopy(smData_hh[hw_s].ProjectionY('scData0_ch%i_h'%(chan), chan+1, chan+1, "e"))
            sampleMean = scData_h.GetMean()
            sampleNoise = scData_h.GetRMS()
            gaus_f = r.TF1('gaus_f','gaus', sampleMean-4*sampleNoise, sampleMean+4*sampleNoise)
            gaus_f.SetParameter(0, 10.0)
            gaus_f.SetParameter(1, sampleMean)
            gaus_f.SetParameter(2, sampleNoise)
            scData_h.Fit(gaus_f, 'QR')
            mean = gaus_f.GetParameter(1)
            rms = gaus_f.GetParameter(2)
            means[chan] = mean
            print mean
            rmss[chan] = rms
            print rms
            pass
        everything[sw_f][0].update(means)
        everything[sw_f][1].update(rmss)
        
        #nsw_f=sw_f
        #nIndices.append(sw_f)

        chList=[]
        meanList=[]
        rmsList=[]
        null=[]
        for chan in xrange(640):
            if layer<5 and chan>512: continue
            chList.append(float(chan))
            meanList.append(float(everything[sw_f][0][chan]))
            rmsList.append(float(everything[sw_f][1][chan]))
            null.append(0)
            pass

        mean_g = r.TGraphErrors( len(chList), np.array(chList), np.array(meanList), np.array(null), np.array(rmsList) )
        #mean_g.SetName('baseline_%i.%i_ge'%(feb,hybrid))
        #mean_g.SetTitle('Baseline vs Channel for %i.%i;Channel;Baseline'%(feb,hybrid))
        mean_g.SetName('baseline_%s_ge'%(hw_s))
        mean_g.SetTitle('Baseline vs Channel for %s;Channel;Baseline'%(hw_s))
        mean_g.SetMarkerStyle(3)
        mean_g_dict[hw_s]=mean_g

        rms_g = r.TGraph( len(chList), np.array(chList), np.array(rmsList) )
        #rms_g.SetName('ENC_%i.%i_g'%(feb,hybrid))
        #rms_g.SetTitle('ENCs vs Channel for %i.%i;Channel;ENC'%(feb,hybrid))
        rms_g.SetName('ENC_%s_g'%(hw_s))
        rms_g.SetTitle('ENCs vs Channel for %s;Channel;ENC'%(hw_s))
        rms_g.SetMarkerStyle(3)
        rms_g_dict[hw_s]=rms_g
        pass
    pass

if options.inRun==0: 
    outFile = r.TFile(options.outfilename,"RECREATE")
else: 
    outFile=r.TFile("fits/hpssvt_%06d_baselineFits.root"%options.inRun, "RECREATE")
outFile.cd()

fDirs={}
nIndices.sort()
print nIndices

fDirs['hh']=outFile.mkdir('hh')
fDirs['baseline']=outFile.mkdir('baseline')
fDirs['enc']=outFile.mkdir('enc')
for feb in xrange(10):
    for hyb in xrange(4):
        hw_s='F%iH%i'%(feb,hyb)
        print "Writing ", hw_s
        fDirs['hh'].cd()
        smData_hh[hw_s].Write()
        fDirs['baseline'].cd()
        mean_g_dict[hw_s].Write()
        fDirs['enc'].cd()
        rms_g_dict[hw_s].Write()

outFile.Close()
inFile.Close()

exit

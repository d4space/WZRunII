#!/bin/bash
#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
#NTUPLEDIR=/data/blue/ksung/TagAndProbeExample
NTUPLEDIR=/u/user/sangilpark/WorkDir/WZ_13TeV/CMSSW_7_2_0/src/WZ/CSA14/TagAndProbe/Efficiency/Zmm_PU4BX50_MuIDISO
#NTUPLEDIR=/u/user/sangilpark/WorkDir/WZ_13TeV/CMSSW_7_2_0/src/WZ/CSA14/TagAndProbe/Efficiency/Zee_PU4BX50_CSA14_EleIDISO

# Z->mumu
#root -l -q plotEff.C+\(\"Muon.bins\",0,0,0,0,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)
root -l -q plotEff.C+\(\"Muon.bins\",1,1,1,1,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)
#root -l -q plotEff.C+\(\"Muon.bins\",1,2,1,2,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)
#root -l -q plotEff.C+\(\"Muon.bins\",1,3,1,3,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)
#root -l -q plotEff.C+\(\"Muon.bins\",1,4,1,4,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)

# Z->ee
#root -l -q plotEff.C+\(\"example.bins\",0,0,0,0,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)
#root -l -q plotEff.C+\(\"example.bins\",1,1,1,1,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)

rm -f *.so *.d

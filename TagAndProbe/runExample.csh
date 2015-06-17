#!/bin/csh

#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

# directory of probes ntuples
set NTUPLEDIR = /u/user/sangilpark/WorkDir/WZ_13TeV/CMSSW_7_2_0/src/WZ/CSA14/TagAndProbe/Efficiency/Zmm_PU4BX50_MuIDISO/
 
# Muon IDISO efficiencies miniAOD sample
#
root -l -q plotEff.C+\(\"example.bins\",0,0,0,0,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)
#root -l -q plotEff.C+\(\"example.bins\",1,1,1,1,\"${NTUPLEDIR}/probes.root\",\"MC\",\"png\",1,0,0\)


rm -f *.so *.d

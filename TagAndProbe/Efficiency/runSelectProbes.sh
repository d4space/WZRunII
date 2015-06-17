#! /bin/bash

#INPUTDIR=/data/blue/ksung/EWKAna/8TeV/Selection
#OUTPUTDIR=/data/blue/ksung/EWKAna/8TeV/Efficiency

INPUTDIR=/d3/scratch/khakim/RunII/WZnTuples/
OUTPUTDIR=/u/user/sangilpark/WorkDir/WZ_13TeV/CMSSW_7_2_0/src/WZ/CSA14/TagAndProbe/Efficiency

#
# Select probes for muon efficiencies
#
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/selectZmm_4bx50.root\",\"${OUTPUTDIR}/Zmm_PU4BX50_MuIDISO\",0,1\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/selectZee_4bx50_CSA14.root\",\"${OUTPUTDIR}/Zee_PU4BX50_CSA14_EleIDISO\",0,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/selectZee_20bx25_CSA14.root\",\"${OUTPUTDIR}/Zee_PU20BX25_CSA14_EleIDISO\",0,1\)

#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff\",0,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuSelEff\",1,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff\",2,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuStaEff\",3,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso\",4,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuPOGIDEff\",5,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuPOGIsoEff\",6,1\)

#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_MuHLTEff\",0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_MuSelEff\",1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_MuTrkEff\",2,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_MuStaEff\",3,0\)

#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuHLTEff\",0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuSelEff\",1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuTrkEff\",2,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuStaEff\",3,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuStaEff_iso\",4,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuPOGIDEff\",5,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_MuPOGIsoEff\",6,0\)

#
# Select probes for electron efficiencies
#
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",0,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",1,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleGsfEff\",2,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",3,1\)

#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_EleHLTEff\",0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_EleSelEff\",1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_pr_select.root\",\"${OUTPUTDIR}/PR_EleGsfEff\",2,0\)

#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_EleHLTEff\",0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_EleSelEff\",1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_EleGsfEff\",2,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_m23_select.root\",\"${OUTPUTDIR}/May23_EleGsfSelEff\",3,0\)

rm *.so *.d

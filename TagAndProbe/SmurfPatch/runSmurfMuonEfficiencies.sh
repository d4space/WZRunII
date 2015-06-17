#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

#RUNLIST=Cert_190456-195658_8TeV_PromptReco_Collisions12_JSON
#TAG=V00-02-04
#MU2012=/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/merged_${RUNLIST}.root
#MCMM=/smurf/dlevans/LeptonTree/V00-02-02/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged_mm.root
#OUTNAME="2012June12V0"

RUNLIST=Cert_190456-196509_8TeV_PromptReco_Collisions12_JSON
TAG=V00-02-06
MU2012=/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/merged_${RUNLIST}.root
MCMM=/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/merged_mm.root
OUTNAME="2012June22V1"
PU=7

#
# options
#

MuonTagAndProbe=1
MuonTagAndProbeMC=5

#
# muons
#

MuonIsoDenominator2012=10
MuonIDDenominator2012=11
MuonFO=12
MuonIDIsoDenominator=13
MuonIDIsoRndTagDenominator=14
MuonIDIsoRndTagDblNoDzDenominator=15

MuonIsoNumerator2012=20
MuonIDNumerator2012=21
MuonIDIsoNumerator2012=22
MuonTrigSglNumerator2012=23
MuonTrigDblLeadNumerator2012=24
MuonTrigDblTrailNumerator2012=25
MuonTrigDblNumerator2012=26

# ID wrt Iso
#root -b -q plotEff.C+\(\"MuonOffline.bins\",2,0,2,1,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_ID\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIsoDenominator2012,$MuonIDNumerator2012\,60,120,\"$MCMM\"\)
#root -b -q plotEff.C+\(\"MuonOffline.bins\",0,0,0,0,\"${MCMM}\",\"HWW_MuonMC${OUTNAME}_NM1Eff_ID\",\"png\",1,${PU},0,$MuonTagAndProbeMC,$MuonIsoDenominator2012,$MuonIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_ID/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_ID/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_ID/eff.root\"\)

# Iso wrt ID
#root -b -q plotEff.C+\(\"MuonOffline.bins\",2,0,2,1,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_Iso\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDDenominator2012,$MuonIsoNumerator2012\,60,120,\"$MCMM\"\)
#root -b -q plotEff.C+\(\"MuonOffline.bins\",0,0,0,0,\"${MCMM}\",\"HWW_MuonMC${OUTNAME}_NM1Eff_Iso\",\"png\",1,${PU},0,$MuonTagAndProbeMC,$MuonIDDenominator2012,$MuonIsoNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_Iso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_MuonMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_Iso/eff.root\"\)

# Single Trigger wrt ID+Iso
root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoDenominator,$MuonTrigSglNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
#root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblLeadNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
#root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblTrailNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)

# Double Trigger (last step = dz cut) wrt ID+Iso + trailing leg (before dz cut) (random tag)
#root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012}\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDblNoDzDenominator,$MuonTrigDblNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"null\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\",\"HWW_Muon${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)

#rm *.so *.d


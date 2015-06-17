#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

#RUNLIST=Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON
#TAG=V00-02-03
#ELE2012=/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/merged_${RUNLIST}.root
#MCEE=/smurf/dlevans/LeptonTree/V00-02-02/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged_ee.root
#OUTNAME="2012May25V0"

#RUNLIST=Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2
#TAG=V00-02-04
#MCEE=/smurf/dlevans/LeptonTree/V00-02-02/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged_ee.root
#ELE2012=/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/merged_${RUNLIST}.root
#OUTNAME="2012June08V1"

#
# v4 with June 12 good run list
#

#RUNLIST=Cert_190456-195658_8TeV_PromptReco_Collisions12_JSON
#TAG=V00-02-04
#MCEE=/smurf/dlevans/LeptonTree/V00-02-02/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_big/merged_ee.root
#ELE2012=/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/merged_${RUNLIST}.root
#OUTNAME="2012June12V0"

RUNLIST=Cert_190456-196509_8TeV_PromptReco_Collisions12_JSON
TAG=V00-02-06
MCEE=/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/merged_ee.root
ELE2012=/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/merged_${RUNLIST}.root
OUTNAME="2012June22V0"
PU=7

#
# options
#

ElectronTagAndProbe=3
ElectronFakeRate=4
ElectronTagAndProbeMC=6

#
# electrons
#

ElectronRecoDenominator=999

ElectronIsoDenominator2012=30
ElectronIDDenominator2012=31
ElectronFO=32
ElectronIDIsoDenominator=33
ElectronIDIsoRndTagDenominator=34
ElectronIDIsoRndTagDblNoDzDenominator=35
ElectronIsoFODenominator2012=36

ElectronIsoNumerator2012=40
ElectronIDNumerator2012=41
ElectronIDIsoNumerator2012=42
ElectronTrigSglNumerator2012=43
ElectronTrigDblLeadNumerator2012=44
ElectronTrigDblTrailNumerator2012=45
ElectronTrigDblNumerator2012=46
ElectronFONumerator=47

#
# alternate factorisation
#

# Isolation wrt FO + ID is as usual

# FO wrt Iso
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,0,2,1,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_FOwrtIso\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIsoDenominator2012,$ElectronFONumerator\,60,120,\"$MCEE\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MCEE}\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_FOwrtIso\",\"png\",1,6,0,$ElectronTagAndProbeMC,$ElectronIsoDenominator2012,$ElectronFONumerator,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_FOwrtIso/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_FOwrtIso/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_FOwrtIso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_FOwrtIso/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_FOwrtIso/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_FOwrtIso/eff.root\"\)

# ID wrt Iso + FO
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,0,2,1,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_IDwrtIsoFO\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIsoFODenominator2012,$ElectronIDNumerator2012\,60,120,\"$MCEE\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MCEE}\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_IDwrtIsoFO\",\"png\",1,6,0,$ElectronTagAndProbeMC,$ElectronIsoFODenominator2012,$ElectronIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_IDwrtIsoFO/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_IDwrtIsoFO/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_IDwrtIsoFO/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_IDwrtIsoFO/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_IDwrtIsoFO/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_IDwrtIsoFO/eff.root\"\)

#
# standard measurements
#

# ID wrt Iso
##root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,0,2,1,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_ID\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIsoDenominator2012,$ElectronIDNumerator2012\,60,120,\"$MCEE\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MCEE}\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_ID\",\"png\",1,${PU},0,$ElectronTagAndProbeMC,$ElectronIsoDenominator2012,$ElectronIDNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_ID/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_ID/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_ID/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_ID/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_ID/eff.root\"\)

# Iso wrt ID
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",2,0,2,1,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_Iso\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDDenominator2012,$ElectronIsoNumerator2012\,60,120,\"$MCEE\"\)
#root -b -q plotEff.C+\(\"ElectronOffline.bins\",0,0,0,0,\"${MCEE}\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_Iso\",\"png\",1,${PU},0,$ElectronTagAndProbeMC,$ElectronIDDenominator2012,$ElectronIsoNumerator2012,60,120\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_Iso/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_Iso/extra\",\"HWW_ElectronMC${OUTNAME}_NM1Eff_Iso/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_Iso/eff.root\"\)

# Single Trigger wrt ID+Iso
root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoDenominator,$ElectronTrigSglNumerator2012\,81,101\)
root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)
root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblLeadNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblTrailNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigTrailDbl/eff.root\"\)

# Double Trigger (last step = dz cut) wrt ID+Iso + trailing leg (before dz cut) (random tag)
#root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012}\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDblNoDzDenominator,$ElectronTrigDblNumerator2012\,81,101\)
#root -b -q plotDataMC.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"null\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)
#root -b -q printEff.C+\(\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/extra\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\",\"HWW_Electron${OUTNAME}_NM1Eff_TrigDzDbl/eff.root\"\)

#rm *.so *.d


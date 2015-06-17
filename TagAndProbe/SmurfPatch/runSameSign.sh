#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

#
# directory of probes ntuples
#

RUNLIST=Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON
TAG=V00-02-07
ELE2012A=/smurf/fgolf/LeptonTree/${TAG}/SingleElectronRun2012/merged_${RUNLIST}.root
MU2012A=/smurf/fgolf/LeptonTree/${TAG}/SingleMuRun2012/merged_${RUNLIST}.root
MCEE=/smurf/fgolf/LeptonTree/V00-02-07/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/merged_ee.root
MCMM=/smurf/fgolf/LeptonTree/V00-02-07/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/merged_mm.root

#
# options
#

PUOPT=5     # PU re-weighting target file number (for MC)
ABSETA=1    # Take the absolute value of eta

MuonTagAndProbe=1
ElectronTagAndProbe=3
MuonTagAndProbeMC=5
ElectronTagAndProbeMC=6

#
# do the measurements for electrons
#
    
ElectronIsoDenominator=30
ElectronIDDenominator=31
ElectronIDIsoDenominator=32
ElectronIDIsoRndTagDenominator=33
ElectronIDIsoRndTagDblNoDzDenominator=34
    
ElectronIsoNumerator=40
ElectronIDNumerator=41
ElectronTrigSglNumerator=42
ElectronTrigDblLeadNumerator=43
ElectronTrigDblTrailNumerator=44
ElectronTrigDblNumerator=45

# ID WRT ISO in EB
 root -b -q plotEff.C+\(\"SameSignElectronEB.bins\",3,1,3,1,\"${ELE2012A}\",\"SameSign_Data2012A_ElectronID_EB\",\"png\",$ABSETA,0,0,$ElectronTagAndProbe,$ElectronIDDenominator,$ElectronIDNumerator\)
 root -b -q plotEff.C+\(\"SameSignElectronEB.bins\",0,0,0,0,\"${MCEE}\",\"SameSign_DY52X_ElectronID_EB\",\"png\",$ABSETA,$PUOPT,0,$ElectronTagAndProbeMC,$ElectronIDDenominator,$ElectronIDNumerator\)
 root -b -q plotDataMC.C+\(\"SameSign_Data2012A_ElectronID_EB/extra\",\"SameSign_DY52X_ElectronID_EB/eff.root\",\"SameSign_Data2012A_ElectronID_EB/eff.root\"\)
 root -b -q printEff.C+\(\"SameSign_Data2012A_ElectronID_EB/extra\",\"SameSign_DY52X_ElectronID_EB/eff.root\",\"SameSign_Muon2012A_ElectronID_EB/eff.root\"\)

# ISO WRT ID in EB
# root -b -q plotEff.C+\(\"SameSignElectronEB.bins\",3,1,3,1,\"${ELE2012A}\",\"SameSign_Data2012A_ElectronIso_EB\",\"png\",$ABSETA,0,0,$ElectronTagAndProbe,$ElectronIsoDenominator,$ElectronIsoNumerator\)
# root -b -q plotEff.C+\(\"SameSignElectronEB.bins\",0,0,0,0,\"${MCEE}\",\"SameSign_DY52X_ElectronIso_EB\",\"png\",$ABSETA,$PUOPT,0,$ElectronTagAndProbeMC,$ElectronIsoDenominator,$ElectronIsoNumerator\)
# root -b -q plotDataMC.C+\(\"SameSign_Data2012A_ElectronIso_EB/extra\",\"SameSign_DY52X_ElectronIso_EB/eff.root\",\"SameSign_Data2012A_ElectronIso_EB/eff.root\"\)
# root -b -q printEff.C+\(\"SameSign_Data2012A_ElectronIso_EB/extra\",\"SameSign_DY52X_ElectronIso_EB/eff.root\",\"SameSign_Muon2012A_ElectronIso_EB/eff.root\"\)

# ID WRT ISO in EE
# root -b -q plotEff.C+\(\"SameSignElectronEE.bins\",3,1,3,1,\"${ELE2012A}\",\"SameSign_Data2012A_ElectronID_EE\",\"png\",$ABSETA,0,0,$ElectronTagAndProbe,$ElectronIDDenominator,$ElectronIDNumerator\)
# root -b -q plotEff.C+\(\"SameSignElectronEE.bins\",0,0,0,0,\"${MCEE}\",\"SameSign_DY52X_ElectronID_EE\",\"png\",$ABSETA,$PUOPT,0,$ElectronTagAndProbeMC,$ElectronIDDenominator,$ElectronIDNumerator\)
# root -b -q plotDataMC.C+\(\"SameSign_Data2012A_ElectronID_EE/extra\",\"SameSign_DY52X_ElectronID_EE/eff.root\",\"SameSign_Data2012A_ElectronID_EE/eff.root\"\)
# root -b -q printEff.C+\(\"SameSign_Data2012A_ElectronID_EE/extra\",\"SameSign_DY52X_ElectronID_EE/eff.root\",\"SameSign_Muon2012A_ElectronID_EE/eff.root\"\)

# ISO WRT ID in EE
# root -b -q plotEff.C+\(\"SameSignElectronEE.bins\",3,1,3,1,\"${ELE2012A}\",\"SameSign_Data2012A_ElectronIso_EE\",\"png\",$ABSETA,0,0,$ElectronTagAndProbe,$ElectronIsoDenominator,$ElectronIsoNumerator\)
# root -b -q plotEff.C+\(\"SameSignElectronEE.bins\",0,0,0,0,\"${MCEE}\",\"SameSign_DY52X_ElectronIso_EE\",\"png\",$ABSETA,$PUOPT,0,$ElectronTagAndProbeMC,$ElectronIsoDenominator,$ElectronIsoNumerator\)
# root -b -q plotDataMC.C+\(\"SameSign_Data2012A_ElectronIso_EE/extra\",\"SameSign_DY52X_ElectronIso_EE/eff.root\",\"SameSign_Data2012A_ElectronIso_EE/eff.root\"\)
# root -b -q printEff.C+\(\"SameSign_Data2012A_ElectronIso_EE/extra\",\"SameSign_DY52X_ElectronIso_EE/eff.root\",\"SameSign_Muon2012A_ElectronIso_EE/eff.root\"\)

# Single Trigger wrt ID+Iso
#echo root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"SameSign_Data2012A_Electron_TrigSgl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoDenominator,$ElectronTrigSglNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Electron_TrigSgl/extra\",\"null\",\"SameSign_Data2012A_Electron_TrigSgl/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Electron_TrigSgl/extra\",\"SameSign_Data2012A_Electron_TrigSgl/eff.root\",\"SameSign_Data2012A_Electron_TrigSgl/eff.root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
#echo root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"SameSign_Data2012A_Electron_TrigLeadDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblLeadNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Electron_TrigLeadDbl/extra\",\"null\",\"SameSign_Data2012A_Electron_TrigLeadDbl/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Electron_TrigLeadDbl/extra\",\"SameSign_Data2012A_Electron_TrigLeadDbl/eff.root\",\"SameSign_Data2012A_Electron_TrigLeadDbl/eff.root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
#echo root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"SameSign_Data2012A_Electron_TrigTrailDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDenominator,$ElectronTrigDblTrailNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Electron_TrigTrailDbl/extra\",\"null\",\"SameSign_Data2012A_Electron_TrigTrailDbl/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Electron_TrigTrailDbl/extra\",\"SameSign_Data2012A_Electron_TrigTrailDbl/eff.root\",\"SameSign_Data2012A_Electron_TrigTrailDbl/eff.root\"\)

# Double Trigger (last step = dz cut) wrt ID+Iso + trailing leg (before dz cut) (random tag)
#echo root -b -q plotEff.C+\(\"ElectronTrigger.bins\",0,0,0,0,\"${ELE2012A}\",\"SameSign_Data2012A_Electron_TrigDzDbl\",\"png\",1,0,0,$ElectronTagAndProbe,$ElectronIDIsoRndTagDblNoDzDenominator,$ElectronTrigDblNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Electron_TrigDzDbl/extra\",\"null\",\"SameSign_Data2012A_Electron_TrigDzDbl/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Electron_TrigDzDbl/extra\",\"SameSign_Data2012A_Electron_TrigDzDbl/eff.root\",\"SameSign_Data2012A_Electron_TrigDzDbl/eff.root\"\)


#
# do the measurements for muons
#

MuonIsoDenominator=10
MuonIDDenominator=11
MuonIDIsoDenominator=12
MuonIDIsoRndTagDenominator=13
MuonIDIsoRndTagDblNoDzDenominator=14
    
MuonIsoNumerator=20
MuonIDNumerator=21
MuonTrigSglNumerator=22
MuonTrigDblLeadNumerator=23
MuonTrigDblTrailNumerator=24
MuonTrigDblNumerator=25

# ID WRT ISO
#echo root -b -q plotEff.C+\(\"SameSignMuon.bins\",4,0,4,1,\"${MU2012A}\",\"SameSign_Data2012A_MuonID\",\"png\",$ABSETA,0,0,$MuonTagAndProbe,$MuonIDDenominator,$MuonIDNumerator,\"${MC}\"\)
#echo root -b -q plotEff.C+\(\"SameSignMuon.bins\",0,0,0,0,\"${MCMM}\",\"SameSign_DY52X_MuonID\",\"png\",$ABSETA,$PUOPT,0,$MuonTagAndProbeMC,$MuonIDDenominator,$MuonIDNumerator\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_MuonID/extra\",\"SameSign_DY52X_MuonID/eff.root\",\"SameSign_Data2012A_MuonID/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Muon_TrigSgl/extra\",\"SameSign_Data2012A_Muon_TrigSgl/eff.root\",\"SameSign_Data2012A_Muon_TrigSgl/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Muon_MuonID/extra\",\"SameSign_DY52X_Muon_MuonID/eff.root\",\"SameSign_Data2012A_Muon_MuonID/eff.root\"\)

# ISO WRT ID
#echo root -b -q plotEff.C+\(\"SameSignMuon.bins\",2,1,2,1,\"${MU2012A}\",\"SameSign_Data2012A_MuonIso\",\"png\",$ABSETA,0,0,$MuonTagAndProbe,$MuonIsoDenominator,$MuonIsoNumerator,\"${MC}\"\)
#echo root -b -q plotEff.C+\(\"SameSignMuon.bins\",0,0,0,0,\"${MCMM}\",\"SameSign_DY52X_MuonIso\",\"png\",$ABSETA,$PUOPT,0,$MuonTagAndProbeMC,$MuonIsoDenominator,$MuonIsoNumerator\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_MuonIso/extra\",\"SameSign_DY52X_MuonIso/eff.root\",\"SameSign_Data2012A_MuonIso/eff.root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_MuonIso/extra\",\"SameSign_DY52X_MuonIso/eff.root\",\"SameSign_Data2012A_MuonIso/eff.root\"\)

# Single Trigger wrt ID+Iso
#echo root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"SameSign_Data2012A_Muon_TrigSgl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoDenominator,$MuonTrigSglNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Muon_TrigSgl/extra\",\"null\",\"SameSign_Data2012A_Muon_TrigSgl/eff.#echo root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Muon_TrigSgl/extra\",\"SameSign_Data2012A_Muon_TrigSgl/eff.#echo root\",\"SameSign_Data2012A_Muon_TrigSgl/eff.#echo root\"\)

# Double Trigger (leading leg) wrt ID+Iso (random tag)
#echo root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"SameSign_Data2012A_Muon_TrigLeadDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblLeadNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Muon_TrigLeadDbl/extra\",\"null\",\"SameSign_Data2012A_Muon_TrigLeadDbl/eff.#echo root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Muon_TrigLeadDbl/extra\",\"SameSign_Data2012A_Muon_TrigLeadDbl/eff.#echo root\",\"SameSign_Data2012A_Muon_TrigLeadDbl/eff.#echo root\"\)

# Double Trigger (trailing leg) wrt ID+Iso (random tag)
#echo root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"SameSign_Data2012A_Muon_TrigTrailDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDenominator,$MuonTrigDblTrailNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Muon_TrigTrailDbl/extra\",\"null\",\"SameSign_Data2012A_Muon_TrigTrailDbl/eff.#echo root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Muon_TrigTrailDbl/extra\",\"SameSign_Data2012A_Muon_TrigTrailDbl/eff.#echo root\",\"SameSign_Data2012A_Muon_TrigTrailDbl/eff.#echo root\"\)

# Double Trigger (last step = dz cut) wrt ID+Iso + trailing leg (before dz cut) (random tag)
#echo root -b -q plotEff.C+\(\"MuonTrigger.bins\",0,0,0,0,\"${MU2012A}\",\"SameSign_Data2012A_Muon_TrigDzDbl\",\"png\",1,0,0,$MuonTagAndProbe,$MuonIDIsoRndTagDblNoDzDenominator,$MuonTrigDblNumerator\,81,101\)
#echo root -b -q plotDataMC.C+\(\"SameSign_Data2012A_Muon_TrigDzDbl/extra\",\"null\",\"SameSign_Data2012A_Muon_TrigDzDbl/eff.#echo root\"\)
#echo root -b -q printEff.C+\(\"SameSign_Data2012A_Muon_TrigDzDbl/extra\",\"SameSign_Data2012A_Muon_TrigDzDbl/eff.#echo root\",\"SameSign_Data2012A_Muon_TrigDzDbl/eff.#echo root\"\)

#rm *.so *.d

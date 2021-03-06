# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerObjects, triggerObjectLabel  = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events("file:/u/user/khakim/WAnalysis/METstudy/CMSSW_7_4_1_patch1/src/WZ/CSA14/test/datasample/DYJetsToLNu_M50_13TeV_bx25sample_MINIAODSIM.root")
events = Events("root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/38D1C54C-0F02-E511-A54E-AC853D9F5256.root")

for iev,event in enumerate(events):
    event.getByLabel(triggerBitLabel, triggerBits)
    event.getByLabel(triggerObjectLabel, triggerObjects)
    event.getByLabel(triggerPrescaleLabel, triggerPrescales)

    print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    print "\n === TRIGGER PATHS ==="
    names = event.object().triggerNames(triggerBits.product())
    for i in xrange(triggerBits.product().size()):
        print "Trigger ", names.triggerName(i), ", prescale ", triggerPrescales.product().getPrescaleForIndex(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)") 

    print "\n === TRIGGER OBJECTS ==="
    for j,to in enumerate(triggerObjects.product()):
        to.unpackPathNames(names);
        print "Trigger object pt %6.2f eta %+5.3f phi %+5.3f  " % (to.pt(),to.eta(),to.phi())
        print "         collection: ", to.collection()
        print "         type ids: ", ", ".join([str(f) for f in to.filterIds()])
        print "         filters: ", ", ".join([str(f) for f in to.filterLabels()])
        pathslast = set(to.pathNames(True))
        print "         paths:   ", ", ".join([("%s*" if f in pathslast else "%s")%f for f in to.filterLabels()]) 
    if iev > 10: break

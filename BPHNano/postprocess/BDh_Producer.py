# This is an example of a NanoAODTools Module to add one variable to nanoAODs.
# Note that:
# -the new variable will be available for use in the subsequent modules
# -it is possible to update the value for existing variables
#
# Example of using from command line:
# nano_postproc.py outDir /eos/cms/store/user/andrey/f.root -I PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule exampleModuleConstr
#
# Example of running in a python script: see test/example_postproc.py
#

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np

class BDhProducer(Module):
    def __init__(self, BSelection):
        self.BSel = BSelection
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("bestB", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        Bcand = Collection(event, "B")
        my_array =np.array(event.B_xgb_2022)
        max_index = np.argmax(my_array)
        self.out.fillBranch("bestB", max_index)
        #eventSum = ROOT.TLorentzVector()
        #for j in filter(self.BSel, Bcand):
        #    eventSum += j.p4()
        #print(len(Bcand))
        #self.out.fillBranch("EventMass", eventSum.M())
        if len(Bcand)>0:
            return True
        else:
            return False


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

BdhModuleConstr = lambda: BDhProducer(BSelection=lambda B: B.pt > 0.5)

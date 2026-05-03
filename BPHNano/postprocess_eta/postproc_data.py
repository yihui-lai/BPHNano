#!/usr/bin/env python3
#
# Example of running the postprocessor to skim events with a cut, and 
# adding a new variable using a Module.
#
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

fnames = ["out_step1.root"]

p = PostProcessor(outputDir=".",
                  inputFiles=fnames,
                  cut="(nEtaTo2L2Pi>=1 || Sum$(EtaMuMu_fitted_mass<0.9 && EtaMuMu_fitted_mass>0.45)>=1 || nEtaTo4Mu>=1)",
                  modules=[],
                  provenance=True,
                  maxEntries=5000000, #just read the first maxEntries events
                  )
p.run()

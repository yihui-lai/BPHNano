#!/usr/bin/env python3
#
# Example of running the postprocessor to skim events with a cut, and 
# adding a new variable using a Module.
#
from BDh_Producer import *

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

fnames = ["test_mc.root"]

p = PostProcessor(outputDir=".",
                  inputFiles=fnames,
                  cut="",
                  modules=[BdhModuleConstr()],
                  provenance=True,
                  maxEntries=50000000, #just read the first maxEntries events
                  )
p.run()

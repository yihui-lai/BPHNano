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

fnames = ["test_data.root"]
#fnames = []
#fnames = [
#        "/eos/cms///store/group/phys_b2g/sqian/VV_comb_workdir/NanoPost/ParkingDoubleMuonLowMass0/BDh_NanoPost_2022_Data_Jun4_ParkingDoubleMuonLowMass0/250611_123452/0000/test_data_651.root",
#        "/eos/cms///store/group/phys_b2g/sqian/VV_comb_workdir/NanoPost/ParkingDoubleMuonLowMass0/BDh_NanoPost_2022_Data_Jun4_ParkingDoubleMuonLowMass0/250611_123452/0000/test_data_652.root",
#        ]
#
#with open("list2.txt", "r") as file:
#    lines = file.readlines()
#    inline = 0
#    for line in lines:
#        inline+=1
#        fnames.append(line.replace("\n",""))
#
#print(fnames)

p = PostProcessor(outputDir=".",
                  inputFiles=fnames,
                  #cut="nLambdabToLambdaMuMu>=1 || nLambdabToLambdahh>=1",
                  #cut="nLambdabToLambdahh>=1",
                  #cut="nEtaTo2L2Pi>=1",
                  cut="nB>=1",
                  modules=[],
                  #modules=[BdhModuleConstr()],
                  provenance=True,
                  maxEntries=5000000, #just read the first maxEntries events
                  )
p.run()

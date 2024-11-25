# nanoAOD producer customized for BPH analysis 

The focus is on B -> J/Psi X analyses.

## Getting started

```shell
cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src
cmsenv
git cms-init
```
Architecture should be el8_amd64_gcc12

## Add the BPHNano package and build everything

```shell
git clone ssh://git@gitlab.cern.ch:7999/btojpsikshort-cpv/reconstruction/BPHNano.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```
or https equivalent

## To run on a test file

```shell
cd PhysicsTools/BPHNano/test/
cmsenv 
cmsRun run_bphNano_cfg.py
```


# nanoAOD producer customized for BPH analysis 

The focus is on B -> mumu X analyses.

## Getting started

```shell
cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src
cmsenv
git cms-init
```
Architecture should be el8 or el9

## Add the BPHNano package and build everything

```shell
git clone this
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


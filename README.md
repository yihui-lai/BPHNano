# nanoAOD producer customized for BParking analysis 

The focus is on RK/K*/phi analyses.

## Getting started

```shell
cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src
cmsenv
git cms-init
```


## Add the BParkingNano package and build everything

```shell
git clone https://github.com/gkaratha/BPHNano ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```

## To run on a test file

```shell
cd PhysicsTools/BPHNano/test/
cmsenv 
cmsRun run_bphNano_cfg.py
```


# nanoAOD producer customized for BPH analysis 

The focus is on Eta->2mu2pi analyses.
Based on the code of offcial BPHnano

## Getting started

Official production, 2022/2023 MC 

```shell
cmssw-el8
cmsrel CMSSW_13_0_23
cd CMSSW_13_0_23/src
cmsenv
git cms-init
```

## Add the BPHNano package and build everything

```shell
git clone git@github.com:yihui-lai/BPHNano.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
git cms-addpkg PhysicsTools/NanoAODTools
cd PhysicsTools/BPHNano/
git checkout yihui
cd ../../
scram b -j 8
```

## To run on a test file

```shell
cd PhysicsTools/BPHNano/postprocess_eta/
cmsenv 
cmsRun run_mc_Run3Summer22.py

```

## To submit jobs to crab 

```shell
source /cvmfs/cms.cern.ch/common/crab-setup.sh
python3  crab_submit_Data_ERA.py DATA.txt
python3  crab_submit_MC_ERA.py   MC.txt
```


The MC config files are generated based on the official command:

```
Run3Summer22NanoAODv12
cmsDriver.py step1  --fileout out_step1.root  --mc --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_mcRun3_2022_realistic_v5 --step NANO --scenario pp --era Run3 --nThreads 2 --python_filename run_mc_Run3Summer22.py --no_exec

Run3Summer22EENanoAODv12
cmsDriver.py step1  --fileout out_step1.root  --mc --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_mcRun3_2022_realistic_postEE_v6 --step NANO --scenario pp --era Run3 --nThreads 2 --python_filename run_mc_Run3Summer22EE.py --no_exec

Run3Summer23NanoAODv12
cmsDriver.py step1  --fileout out_step1.root  --mc --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_mcRun3_2023_realistic_v15 --step NANO --scenario pp --era Run3_2023 --nThreads 2 --python_filename run_mc_Run3Summer23.py --no_exec

Run3Summer23BPixNanoAODv12
cmsDriver.py step1  --fileout out_step1.root  --mc --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_mcRun3_2023_realistic_postBPix_v6 --step NANO --scenario pp --era Run3_2023 --nThreads 2 --python_filename run_mc_Run3Summer23BPix.py --no_exec

```

The Data config files are generated based on the official command:

```
Run3Summer22NanoAODv12 (CDE)
cmsDriver.py step1  --fileout out_step1.root  --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_dataRun3_v2 --step NANO --scenario pp --era Run3 --nThreads 2 --python_filename run_data_Run3Summer22.py --no_exec --data

Run3Summer22EENanoAODv12 (FG)
cmsDriver.py step1  --fileout out_step1.root  --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_dataRun3_PromptAnalysis_v1 --step NANO --scenario pp --era Run3 --nThreads 2 --python_filename run_data_Run3Summer22EE.py --no_exec --data

Run3Summer23NanoAODv12
cmsDriver.py step1  --fileout out_step1.root  --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_dataRun3_PromptAnalysis_v1 --step NANO --scenario pp --era Run3_2023 --nThreads 2 --python_filename run_data_Run3Summer23.py --no_exec --data

Run3Summer23BPixNanoAODv12
cmsDriver.py step1  --fileout out_step1.root  --eventcontent NANOAODSIM --datatier NANOAOD --conditions 130X_dataRun3_PromptAnalysis_v1 --step NANO --scenario pp --era Run3_2023 --nThreads 2 --python_filename run_data_Run3Summer23BPix.py --no_exec --data
```




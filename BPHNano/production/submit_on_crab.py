#!/usr/bin/env python

import os
from argparse import ArgumentParser
from fnmatch import fnmatch
import yaml

import re
import datetime

from schema import Schema, And, Or, Optional, SchemaError

import CRABClient

from CRABAPI.RawCommand import crabCommand

from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException 

from CRABClient.UserUtilities import config
from multiprocessing import Process


production_tag = datetime.date.today().strftime('%Y%b%d')


def parse_args():
    parser = ArgumentParser(description="A multicrab submission script")
    parser.add_argument('-y', '--yaml', default = 'test_samples.yml', help = 'File with dataset descriptions')
    parser.add_argument('-c', '--cmd', default='submit', choices = ['submit', 'status'], help= 'Crab command')
    parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed') 
    parser.add_argument('-w', '--workarea', default='BPHNANO_%s' % production_tag, help = 'Crab working area name')
    parser.add_argument('-o', '--outputdir', default= '/store/user/valukash/', help='LFN Output high-level directory: the LFN will be saved in outputdir+workarea ')
    parser.add_argument('-s', '--site', default='T2_CH_CSCS', help='T2 or T3 cite where user has access. To be checked with crab checkout')
    parser.add_argument('-t', '--tag', default=production_tag, help='Production Tag extra')
    parser.add_argument('-p', '--psetcfg', default="../test/run_bphNano_cfg.py", help='Plugin configuration file')
    parser.add_argument('-e', '--extra', nargs='*', default=list(),  help='Optional extra input files')
    parser.add_argument('-tt', '--test', action='store_true', help='Flag a test job')
    return parser.parse_args()
    
def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

def status(directory):
    try:
        crabCommand('status', dir=directory)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))


expected_schema = Schema({
    "common": {
        "data": {
            "lumimask": And(str, error="lumimask should be a string"),
            "splitting": And(int, error="splitting should be an integer"),
            "globaltag": And(str, error="globaltag should be a string"),
        },
        "mc": {
            "splitting": And(int, error="splitting should be an integer"),
            "globaltag": And(str, error="globaltag should be a string"),
        },
    },
    "samples": And(dict, error="samples should be a dict with keys dataset (str), isMC (bool). Optional keys: globaltag (str), parts (list(int))")
    }
    )

samples_schema = Schema({
    "dataset": And(str, error="dataset should be a string"),
    "isMC": And(bool, error="isMC should be a boolean"),
    Optional("decay") : And(str, error="decay to reconstruct"), 
    Optional("goldenjson") : And(str, error="golden json file path should be a string"),
    Optional("globaltag") : And(str, error="globaltag should be a string"),
    Optional("parts"): [And(int, error="parts should be a list of integers")]
})


def validate_yaml(data):
    try:
       expected_schema.validate(data)
       for name, content in data["samples"].items():
           samples_schema.validate(content)
       print("YAML structure is valid.")
    except SchemaError as e:
       print("YAML structure is invalid:", e)
       import sys
       sys.exit(1)

def get_common_config(args):

    config_ = config()
    
    config_.General.transferOutputs = True
    config_.General.transferLogs = True
    config_.General.workArea = args.workarea

    config_.Data.publication = False
    config_.Data.outLFNDirBase = args.outputdir + '/'+ config_.General.workArea
    config_.Data.inputDBS = 'global'

    config_.JobType.pluginName = 'Analysis'
    config_.JobType.psetName = args.psetcfg
    config_.JobType.maxJobRuntimeMin = 2700 #can not use with Automatic 
    config_.JobType.allowUndistributedCMSSW = True
    config_.JobType.inputFiles = args.extra

    config_.Site.storageSite = args.site

    return config_

def get_dataset_config(config, common_config, dataset_config, production_tag):
     isMC = dataset_config['isMC']       
     data_type = 'mc' if isMC else 'data'

     config.Data.splitting = 'FileBased'
     if not isMC:
         config.Data.lumiMask = dataset_config.get('lumimask', None)
     else:
         config.Data.lumiMask = ''
     config.Data.unitsPerJob = common_config[data_type].get('splitting', None)

     globaltag = dataset_config.get('globaltag', "auto:run3_data")
     if globaltag == "auto:run3_data":
         globaltag = common_config[data_type].get('globaltag', "auto:run3_data")

     decay = dataset_config.get('decay', 'KshortLL')
     
     maxevents = -1
        
     config.JobType.pyCfgParams = [
                    'isMC=%s' % isMC, 'reportEvery=1000',
                    'tag=%s' % production_tag,
                    'globalTag=%s' % globaltag,
                    'decay=%s' % decay,
                    'maxEvents=%s' % maxevents,
     ]



     return config
   

if __name__ == '__main__':

    args = parse_args()
    with open(args.yaml, "r") as f:
        samples = yaml.safe_load(f) # Parse YAML file
    validate_yaml(samples)
  
    if args.cmd == "submit":
        print("")
        print(f"Submit Crab jobs for {args.yaml} with filter {args.filter} applied")
        config = get_common_config(args)
        common = samples['common'] if 'common' in samples else {'data' : {}, 'mc' : {}}
        # loop over samples
        for sample, info in samples['samples'].items():
            # Given we have repeated datasets check for different parts
            
            parts = info['parts'] if 'parts' in info else [None]
            for part in parts:
                name = sample % part if part is not None else sample
                config.Data.outLFNDirBase = args.outputdir + config.General.workArea
                config.General.workArea = "workarea/" + args.workarea + "_" + name
                
                # filter names according to what we need
                if not fnmatch(name, args.filter): continue
        
                config.Data.inputDataset = info['dataset'] % part \
                                         if part is not None else \
                                                  info['dataset']
                get_dataset_config(config, common, info, production_tag)
                if args.test:
                   config.Data.totalUnits = 10

                config.General.requestName = name + "_" + production_tag
                config.JobType.outputFiles = ['_'.join(['BPHNano', 'mc' if info['isMC'] else 'data', production_tag])+'.root']
 

                print(f"Submit Crab job for {name}")
                print(config)   
                submit(config)
    elif args.cmd == "status":
        print(f"Getting crab status for {args.dir}")
        status(args.dir)
    else:
        print(f"Invalid Crab command : {args.cmd}")
    


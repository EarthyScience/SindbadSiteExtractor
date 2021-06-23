#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:57:14 2021

@author: simon
"""
#%% Load library
import argparse
import json
import xarray as xr
import multiprocessing as mp
import sys
sys.path.append('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/cliff_workflow/')

#%% Retrieve arguments to parse in data processing routine
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file',
                        '-config_file',
                        type=str,
                        help='Path to the config file',
                        default='/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/cliff_workflow/exp_config.json',
                        required=False)   
    args = parser.parse_args()

    return args

#%% Retrieve parsed arguments
args = parse_args()

#%% Open config files
with open(args.config_file, 'r') as myfile:
    exp_config=myfile.read()
exp_config = json.loads(exp_config)

#%% Export variable information json file
with open(exp_config['OutPath']['nc_file'] + "/hourly/" + exp_config['FLUXNET_version'] + "/variable_info.json", 'w') as fp:
    json.dump(exp_config['variables'], fp, indent=2, sort_keys=True)

#%% Check sites to process
if exp_config["sites_toRun"] == ['all']:
    sites = xr.open_zarr(exp_config["DataPath"]['fluxcom_cube']).site.values
else:
    sites = exp_config["sites_toRun"]
    
#%% Run model run
from utils.cliff_workflow import DataProcessor
print('-------------------------------------------\n'
      'Initialization of the data processing\n'
      '-------------------------------------------')
if __name__ == '__main__':
    exp_settings = {'exp_settings' : exp_config}
    data_process_run = DataProcessor(**exp_settings)
    p = mp.Pool(exp_settings['exp_settings']['njobs_parallel'], maxtasksperchild=1)
    p.map(data_process_run.compile_data, sites)
    p.close()
    p.join()
    print('----------------------------\n'
          'End of model data processing\n'
        '----------------------------')
    

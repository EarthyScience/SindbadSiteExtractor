#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:57:14 2021

@author: sujan
"""
from multiprocessing import set_start_method, Pool
# set_start_method("spawn")
from multiprocessing import get_context


import sys
import os
import utils.shared_utils as shut
import utils.setup_config as set_conf
sys.path.append(os.getcwd())

def get_site_list(cubepath, version):
    import numpy as np
    site_list_path = os.path.join(exp_config['output_dir_path']['info_config'], f"site_list_{version}.csv")
    print(site_list_path)
    if os.path.exists(site_list_path):
        site_list = [site.strip() for site in open(site_list_path).readlines()]
    else:
        from fluxcom.providers import eddy_covariance as ec
        from fluxcom.core.variables import Variable
        Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = cubepath,version  = version)
        Climate_FLUXNET_data = Climate_FLUXNET_prov.get_data(Variable('P'))
        site_list_d = Climate_FLUXNET_data.site.values
        site_list = sorted(site_list_d)
        print(site_list_path)
        np.savetxt(site_list_path, site_list, delimiter=",", fmt="%s")
    return site_list

class DataCompiler:    
    def __init__(self, exp_settings=None):
        self.exp_settings = exp_settings
    def process(self, site):
        site_info = compile_site_data(site, self.exp_settings)
        if isinstance(site_info, dict):
            shut.save_json(os.path.join(exp_config['output_dir_path']['site_info'], f'{site}_info_{exp_config["FLUXNET_version"]}.json'), site_info)
        return 'done'


if __name__ == '__main__':
    #%% Retrieve parsed arguments
    if len(sys.argv) > 1:
        exp_config_path = sys.argv[1]
        if not os.path.isabs(exp_config_path):
            sys.exit(f'The provided path for exp config: {exp_config_path} is not absolute. Provide full absolute path to the configuration file.')
        if len(sys.argv) > 2:
            fn_versions =[sys.argv[2]]
        else:
            fn_versions=[]
        if len(sys.argv) > 3:
            tm_scales =[sys.argv[3]]
        else:
            tm_scales=[]
        if len(sys.argv) > 4:
            output_dir_path =sys.argv[4]
        else:
            output_dir_path=""
    else:
        sys.exit(f'run_extraction cannot extract data without the path of the main extraction config as the only argument. Provide absolute path to the configuration file.')
    # get basic config to check the sources and loop
    exp_config_r = set_conf.get_inp_config(exp_config_path)

    if len(fn_versions) == 0:
        fn_versions = exp_config_r["FLUXNET_version"]
        if isinstance(fn_versions, list):
            fn_versions=[_fn.strip() for _fn in fn_versions if len(_fn) != 0]
        else:
            fn_tmp = fn_versions.split(" ")
            fn_versions=[_fn.strip() for _fn in fn_tmp if len(_fn) != 0]

    if len(tm_scales) == 0:
        tm_scales = exp_config_r["temporal_resolution"]
        if isinstance(tm_scales, list):
            tm_scales=[_tm.strip() for _tm in tm_scales if len(_tm) != 0]
        else:
            tm_tmp = tm_scales.split(" ")
            tm_scales=[_tm.strip() for _tm in tm_tmp if len(_tm) != 0]
    if len(output_dir_path) == 0:
        output_dir_path = exp_config_r["output_dir_path"]

    # loop through sites and get full configurations
    for tm_scale in tm_scales:
        for fn_version in fn_versions:
            exp_config = set_conf.get_exp_configuration(exp_config_path, fn_version=fn_version, temporal_resolution=tm_scale, output_dir_path=output_dir_path)

            # #%% Check sites to process
            if f'sites_{fn_version}' not in locals():
                if exp_config["sites_toRun"] == ['all']:
                    sites = get_site_list(exp_config["fluxcom_cube_path"], exp_config["FLUXNET_version"])
                    # eval(f'sites_{fn_version} = sites')
                else:
                    sites = exp_config["sites_toRun"]
            else:
                sites = eval(f'sites_{fn_version}')

            print(f'{exp_config["FLUXNET_version"]}: {sites}')

            exp_settings = {'exp_settings' : exp_config}
            from utils.workflow_extraction import compile_site_data
            data_process_run = DataCompiler(**exp_settings)
            if  exp_config["njobs_parallel"] > 1:
                p = Pool(exp_settings['exp_settings']['njobs_parallel'], maxtasksperchild=1)
                site_info_dict = p.map(data_process_run.process, sites)
                p.close()
                p.join()
            else:
                site_info_dict ={}
                for site in sites:
                    site_info_dict[site] = compile_site_data(site, **exp_settings)
            print('The extraction is complete. Back to run_extract.py.')

            print('----------------------------\n'
                'End of model data processing\n'
                '----------------------------')

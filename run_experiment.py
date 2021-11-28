#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:57:14 2021

@author: sujan
"""
import json
import sys
import os
import copy
import utils.shared_utils as shut
import utils.setup_config as set_conf
sys.path.append(os.getcwd())

def get_site_list(cubepath, version):
    import numpy as np
    site_list_path = os.path.join(exp_config['OutPath']['info_config'], f"site_list_{version}.csv")
    print(site_list_path)
    if os.path.exists(site_list_path):
        site_list = [site.strip() for site in open(site_list_path).readlines()]
    else:
        from fluxcom.providers import eddy_covariance as ec
        from fluxcom.core.variables import Variable
        Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = cubepath,version  = version)
        Climate_FLUXNET_data = Climate_FLUXNET_prov.get_data(Variable('P'))
        site_list = Climate_FLUXNET_data.site.values
        print(site_list_path)
        np.savetxt(site_list_path, site_list, delimiter=",", fmt="%s")
    return site_list

class DataCompiler:    
    def __init__(self, exp_settings=None):
        self.exp_settings = exp_settings
    def process(self, site):
        site_data = compile_site_data(site, self.exp_settings)
        return site_data


if __name__ == '__main__':
    #%% Retrieve parsed arguments

    if len(sys.argv) > 1:
        config_filepath = sys.argv[1]
    else:
        config_filepath = os.path.join(os.getcwd(),'configuration', 'exp_config.json')

    # get basic config to check the sources and loop
    exp_config_r = set_conf.get_inp_config(config_filepath)

    fn_versions = exp_config_r["FLUXNET_version"]
    tm_scales = exp_config_r["temporal_resolution"]

    if isinstance(fn_versions, list):
        fn_versions=[_fn.strip() for _fn in fn_versions if len(_fn) != 0]
    else:
        fn_tmp = fn_versions.split(" ")
        fn_versions=[_fn.strip() for _fn in fn_tmp if len(_fn) != 0]

    if isinstance(tm_scales, list):
        tm_scales=[_tm.strip() for _tm in tm_scales if len(_tm) != 0]
    else:
        tm_tmp = tm_scales.split(" ")
        tm_scales=[_tm.strip() for _tm in tm_tmp if len(_tm) != 0]

    # loop through sites and get full configurations
    for tm_scale in tm_scales:
        for fn_version in fn_versions:
            exp_config = set_conf.get_exp_configuration(config_filepath, fn_version=fn_version, temporal_resolution=tm_scale)

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
            from utils.workflow_experiment import compile_site_data
            data_process_run = DataCompiler(**exp_settings)
            if  exp_config["njobs_parallel"] > 1:
                import multiprocessing as mp
                p = mp.Pool(exp_settings['exp_settings']['njobs_parallel'], maxtasksperchild=1)
                site_info_dict = p.map(data_process_run.process, sites)
                p.close()
                p.join()
            else:
                site_info_dict ={}
                for site in sites:
                    site_info_dict[site] = compile_site_data(site, **exp_settings)

            if isinstance(site_info_dict, dict):
                shut.save_json(os.path.join(exp_config['OutPath']['info_config'], 'site_info.json'), site_info_dict)
            else:
                shut.save_json(os.path.join(exp_config['OutPath']['info_config'], 'site_info.json'), site_info_dict)

            print('----------------------------\n'
                'End of model data processing\n'
                '----------------------------')

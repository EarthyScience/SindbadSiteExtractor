#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 10:50:58 2021

@author: simon
"""
#%% Load library
from fluxcom.providers import eddy_covariance as ec
from fluxcom.variable import Variable
from fluxcom.providers.transformers import hourly_to_daily
import numpy as np
import pandas as pd
from datetime import datetime
import xarray as xr
import argparse
import json
import glob
import os
import matplotlib.pyplot as plt

#%% Retrieve arguments to parse in data processing routine
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file',
                        '-config_file',
                        type=str,
                        help='Path to the config file',
                        default='/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/output_cliff/exp_config.json',
                        required=False)   
    args = parser.parse_args()

    return args

#%% Retrieve parsed arguments
args = parse_args()

#%% Open config files
with open(args.config_file, 'r') as myfile:
    exp_config=myfile.read()
exp_config = json.loads(exp_config)
flux_source = ['fluxnetBGI2021.BRK15.DD', 'fluxnetBGI2021.LTL07.DD']
clim_source = ['CRUNCEP', 'ERAinterim']

for flux_ in flux_source:

    for clim_ in clim_source:

        #%% Export /fluxnetBGI2021.BRK15.DD/CRUNCEP/ files
        csv_files = glob.glob(exp_config['DataPath']['cliff_output'] + '/' + flux_ + '/' + clim_ + '/*.csv')
        for site_ in csv_files:
            cliff_data = pd.read_csv(site_)
            date_ = np.array([np.datetime64(datetime.strptime(str(np.array(cliff_data['year'])[obs_]) + " " + str(np.array(cliff_data['julday'])[obs_]), '%Y %j').strftime("%Y-%m-%d")) for obs_ in np.arange(cliff_data.shape[0])])
            
            #%% create downscaled xarray
            downscaled_data = []
            clim_vars= list(exp_config['variables'].keys())
            for var_ in clim_vars:
                dat_ = np.array(cliff_data[var_])
                dat_ = xr.DataArray(dat_,dims=['time'], coords={'time': date_}).to_dataset(name=var_)
                unit_scalar = int(exp_config['variables'][var_]['source2sindbadUnit'])
                dat_[var_] = dat_[var_] * unit_scalar
                dat_[var_] = dat_[var_].assign_attrs(units = exp_config['variables'][var_]['variableUnit'],
                                               short_name = exp_config['variables'][var_]['nameShort'],
                                               long_name = exp_config['variables'][var_]['nameLong'])
                downscaled_data.append(dat_)
            downscaled_data = xr.merge(downscaled_data)
            
            #%% Export files
            if not os.path.exists(exp_config['OutPath']['nc_file'] + '/' + flux_ + '/' + clim_):
                os.makedirs(exp_config['OutPath']['nc_file'] + '/' + flux_ + '/' + clim_) 
            downscaled_data.to_netcdf(exp_config['OutPath']['nc_file'] + '/' + flux_ + '/' + clim_ + '/' + os.path.basename(site_).split('.csv')[0] + '.nc', mode= 'w')
            
            #%% Load original data
            if flux_.split('.')[1] == 'BRK15':
                version_ = 'FLUXNET2015'
            elif flux_.split('.')[1] == 'LTL07':
                version_ = 'LaThuile'
            Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath   = "/Net/Groups/BGI/scratch/FLUXCOM/v0.01/site_cube",
                                                                             transforms = hourly_to_daily,
                                                                             version    = version_, 
                                                                             site       = os.path.basename(site_).split('.')[0])
            Climate_FLUXNET_vars = []
            for var_ in clim_vars:
               Climate_FLUXNET_vars.append(Variable(var_))   
            Climate_FLUXNET_data = Climate_FLUXNET_prov.get_data(Climate_FLUXNET_vars)
            Climate_FLUXNET_data = Climate_FLUXNET_data.rename(**{v:v+'_original' for v in list(set(Climate_FLUXNET_data.variables)-set(Climate_FLUXNET_data.coords))})
            
            #%% Define plot matrix
            plot_matrix = {"row_1": {"col_1":{"x_data": "TA_original", "x_axis": "TA original", "y_data": "TA", "y_axis": "TA cliff"},
                                     "col_2":{"x_data": 'P_original', "x_axis": 'P original', "y_data": 'P', "y_axis": 'P cliff'},
                                     "col_3":{"x_data": 'VPD_original', "x_axis": 'VPD original', "y_data": 'VPD', "y_axis": 'VPD cliff'}},                 
                           "row_2":{"col_1":{"x_data": "SW_IN_original", "x_axis": "SW_IN original", "y_data": "SW_IN", "y_axis": "SW_IN cliff"},
                                    "col_2":{"x_data": 'NETRAD_original', "x_axis": 'NETRAD original', "y_data": 'NETRAD', "y_axis": 'NETRAD cliff'}}}

            
            #%% Plot comparison
            if not os.path.exists(exp_config['OutPath']['figs'] + '/' + flux_ + '/' + clim_):
                os.makedirs(exp_config['OutPath']['figs'] + '/' + flux_ + '/' + clim_) 
            plot_data = xr.merge([Climate_FLUXNET_data, downscaled_data])
            fig, ax = plt.subplots(2,3, figsize=(10,6), constrained_layout=True)
            for row_ in range(len(list(plot_matrix))):
                row_class = plot_matrix[list(plot_matrix)[row_]]
                for col_ in range(len(list(row_class))):
                    x_data = plot_data[row_class[list(row_class)[col_]]['x_data']]
                    y_data = plot_data[row_class[list(row_class)[col_]]['y_data']]
                    QA_data = Climate_FLUXNET_prov.get_data(Variable(row_class[list(row_class)[col_]]['y_data'] + '_QC'))
                    if row_class[list(row_class)[col_]]['y_data'] == 'P':
                        mask_QA = (QA_data >=.85) & np.isfinite(x_data) & (x_data>= 0)
                    else:
                        mask_QA = (QA_data >=.85) & (np.isfinite(x_data))
                    ax[row_,col_].scatter(x_data.where(mask_QA), y_data.where(mask_QA), label = 'original data')
                    ax[row_,col_].scatter(x_data.where(~mask_QA), y_data.where(~mask_QA), label = ' bad quality')
                    ax[row_,col_].spines['top'].set_visible(False)
                    ax[row_,col_].spines['right'].set_visible(False)
                    ax[row_,col_].set_xlabel(row_class[list(row_class)[col_]]['x_axis'], size=14)
                    ax[row_,col_].set_ylabel(row_class[list(row_class)[col_]]['y_axis'], size=14)
                    ax[row_,col_].plot(x_data.where(mask_QA), x_data.where(mask_QA), color ='black', 
                                       linestyle='dashed', linewidth=2)
                    max_abs_diff = np.round(np.nanmax(np.abs(x_data.where(mask_QA) -  y_data.where(mask_QA))), 15)
                    ax[row_,col_].set_title(f'max abs diff = {max_abs_diff}')
                    ax[row_,col_].legend(loc='lower right', frameon=False)
            fig.delaxes(ax[1,2])
            plt.savefig(exp_config['OutPath']['figs'] + '/' + flux_ + '/' + clim_ + '/' + os.path.basename(site_).split('.csv')[0] + '.png', dpi=300, transparent=False)
        
            
            
            
            
    
    




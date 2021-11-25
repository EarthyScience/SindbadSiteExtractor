#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 22 2021

@author: sujan
"""
#%% Load libray
import xarray as xr
import numpy as np
import json
import logging
logger = logging.getLogger(__name__)
#%% Process data
def data_structure_temporal_NoDepth(data, lat, lon):

    var_ = list(data.keys())
    # print(data,data.to_array().values.shape)
    ds_ = xr.DataArray(data.to_array().values.reshape(len(var_),-1,1,1), 
                 dims=['variable', 'time', 'latitude', 'longitude'],
                 coords={'variable': var_,
                         'time': data.time.values,
                         'latitude': [lat],
                         'longitude': [lon]}).to_dataset(dim="variable")
    return ds_

def flatten_hour_to_time(_data):
    import pandas as pd
    # print(_data)
    time1=_data['time'].values[0]
    _data = _data.rename({'time':'date'}).stack(time=('date','hour')).reset_index('time')
    times = pd.date_range(time1, periods=len(_data['time']), freq='H')
    _data['time']=times
    return _data

def log_site_info(dataset, site, t_reso, src_var_name, tar_var_name, data_units, json_units_tar, json_units_src, bounds_in, data, transform, partitioning='None'):
    tname = get_transform_name(transform)
    ## get data information
    data_a = data.values
    data_min = np.nanmin(data_a)
    data_max = np.nanmax(data_a)
    ndata = data_a.size
    nan_c = np.isnan(data_a).sum()
    nan_c_p = round(nan_c * 100/(ndata), 2)
    if nan_c_p > 90:
        logger.warning(f'nanwarn: {nan_c_p} % of data points are nan.')
    data_tmp = data_a[~np.isnan(data_a)]
    n_data_tmp = data_tmp.size
    n_out_bounds = ((bounds_in[0] > data_tmp) & (data_tmp > bounds_in[1])).sum()
    n_out_bounds_perc= round(n_out_bounds * 100/(n_data_tmp), 2)
    logger.info(f'\n Dataset:: {dataset} \n site: {site}\
                \n reso:: {t_reso}\
                \n var:: \n \t src: {src_var_name} \n \t target: {tar_var_name}\n \t partitioning: {partitioning} \n units:: \n \t data: {data_units} \n \t json-target:: {json_units_tar} \n \t json-source: {json_units_src}, \n data (after processing):: \n \t shape: {data.values.shape} \n \t size: {data.values.size} \n \t min: {data_min} \n \t max: {data_max} \n json_bnds:: \n \t min: {bounds_in[0]} \n \t max: {bounds_in[1]} \n counts:: \n \t out_bounds_none_NaN: {n_out_bounds} ({n_out_bounds_perc}%) \n \t nans: {nan_c} ({nan_c_p}%) \n transform:: {tname}')
    return


def save_json(json_path, data_to_save):
    with open(json_path, 'w') as fp:
        json.dump(data_to_save, fp, indent=2, sort_keys=True)
    return


def set_units(data, var_name, src_units, tar_units, units_scalar):
    data_out = data.copy()
    if units_scalar == 1.:
        if src_units.lower() != tar_units.lower():
            logger.warning(f'Unit Conversion::\n{var_name}: \nunits \n\t src: [{src_units}] \n\t target: [{tar_units}] \n scaled by: {units_scalar} \n Note :: potentially inconsistent units. Make sure that transforms below take care of unit conversion.')
        else:
            logger.warning(f'Unit Conversion::\n{var_name}: \nunits \n\t src: [{src_units}] \n\t target: [{tar_units}] \n scaled by: {units_scalar} \n Note :: No unit conversion needed.')

            data_out.attrs["units"]=src_units
    else:
        data_out.values = data_out.values * units_scalar
        logger.info(f'Unit Conversion::\n{var_name}: \nunits \n\t src: [{src_units}] \n\t target: [{tar_units}] \n scaled by: {units_scalar} \n mean \n\t before : {np.nanmean(data.values)} \n\t after : {np.nanmean(data_out.values)} ')
        data_out.attrs["units"]=tar_units
    return data_out

def get_transform_name(_transform):
    if _transform is None or _transform is False or len(_transform) == 0:
        tname="None"
    else:
        tname = ",".join([__tr.class_name for __tr in _transform])
    return tname

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Utilities: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('shared functions: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
    
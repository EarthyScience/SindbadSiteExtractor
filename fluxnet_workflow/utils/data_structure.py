#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 10:19:47 2021

@author: simon
"""
#%% Load libray
import xarray as xr

#%% Process data
def data_structure_temporal_NoDepth(data, lat, lon):

    var_ = list(data.keys())
    ds_ = xr.DataArray(data.to_array().values.reshape(len(var_),-1,1,1), 
                 dims=['variable', 'time', 'lat', 'lon'],
                 coords={'variable': var_,
                         'time': data.time.values,
                         'lat': [lat],
                         'lon': [lon]}).to_dataset(dim="variable")
    return ds_


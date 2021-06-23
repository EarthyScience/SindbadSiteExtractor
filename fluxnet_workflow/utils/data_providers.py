#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 10:21:27 2021

@author: simon
"""
from fluxcom.providers import eddy_covariance as ec
from fluxcom.providers.meteo_downscale import era5
from fluxcom.providers.modis.refl import modis_MCD43A
from fluxcom.providers.modis.LST import modis_MxD11A
from fluxcom.providers.transformers import hourly_to_daily
from fluxcom.providers import transformers
from fluxcom.variable import Variable
import xarray as xr
import pandas as pd
import numpy as np
import sys
sys.path.append('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/')
from utils.data_structure import data_structure_temporal_NoDepth

class EddyCovariance_FLUXNET:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.version = config["FLUXNET_version"]
        self.vars = config["variables"]["EddyCovariance_vars"]["FLUXNET_EddyCovariance"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]        
     
    def EddyCovariance_FLUXNET_proc(self):
        
        #%% Process data
        ec_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                                     transforms = hourly_to_daily,
                                                     NEE_partitioning_method = None,
                                                     version    = self.version, 
                                                     site       =self.site)
        ec_vars = []
        for k, v in self.vars.items():
            ec_vars.append(Variable(v['sourceVariableName'], unit = v['sourceVariableUnit'], partitioning = v['partitioning']))   
        ec_data = ec_prov.get_data(ec_vars)
        ec_data = ec_data.rename(**{v:v+'_FLUXNET' for v in list(set(ec_data.variables)-set(ec_data.coords))})
                
        #%% Retrieve weighted uncertainty - NEE
        NEE_unc_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                                          NEE_partitioning_method = None,
                                                          version    = self.version, 
                                                          site       =self.site)
        weighted_NEE_unc_transformer = [transformers.CarbonUncertaintyWeightedHourToDay(NEE_unc_prov.get_data(Variable('NEE_QC')))]
        NEE_unc_prov.add_transform(weighted_NEE_unc_transformer)
        weighted_NEE_unc = NEE_unc_prov.get_data(Variable('NEE_RANDUNC'))   
        weighted_NEE_unc = weighted_NEE_unc.to_dataset(name=weighted_NEE_unc.variable_name).rename({'NEE_RANDUNC': 'NEE_RANDUNC_weighted_FLUXNET'})                                                                         
        
        #%% Retrieve weighted uncertainty - LE
        LE_unc_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                                         NEE_partitioning_method = None,
                                                         version    = self.version, 
                                                         site       =self.site)
        weighted_LE_unc_transformer = [transformers.UncertaintyWeightedHourToDay(LE_unc_prov.get_data(Variable('LE_QC')))]
        LE_unc_prov.add_transform(weighted_LE_unc_transformer)
        weighted_LE_unc = LE_unc_prov.get_data(Variable('LE_RANDUNC'))   
        weighted_LE_unc = weighted_LE_unc.to_dataset(name=weighted_LE_unc.variable_name).rename({'LE_RANDUNC': 'LE_RANDUNC_weighted_FLUXNET'})                                                                         
        
        #%% Retrieve weighted uncertainty - H
        H_unc_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                                        NEE_partitioning_method = None,
                                                        version    = self.version, 
                                                        site       =self.site)
        weighted_H_unc_transformer = [transformers.UncertaintyWeightedHourToDay(H_unc_prov.get_data(Variable('H_QC')))]
        H_unc_prov.add_transform(weighted_H_unc_transformer)
        weighted_H_unc = H_unc_prov.get_data(Variable('H_RANDUNC'))   
        weighted_H_unc = weighted_H_unc.to_dataset(name=weighted_H_unc.variable_name).rename({'H_RANDUNC': 'H_RANDUNC_weighted_FLUXNET'})                                                                         
        ec_data = xr.merge([ec_data, weighted_NEE_unc, weighted_LE_unc, weighted_H_unc])
        
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            ec_data[var_] = ec_data[var_] * unit_scalar
            #if (np.nanmin(ec_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(ec_data[var_]) > self.vars[var_]['bounds'][1]):
            #    print("Values are out of bounds:\n"+var_)
        
        #%% Restructure data
        ec_data = data_structure_temporal_NoDepth(ec_data, self.lat, self.lon)
        
        return ec_data#.sel(time = slice(self.start_date, self.end_date))

class SWC_TS_FLUXNET:
    def __init__(self, site_info, config):
    
        #%% Process data
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.version = config["FLUXNET_version"]
        self.vars = config["variables"]["SWC_TS_vars"]['FLUXNET_SWC_TS']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]   
    
    def SWC_TS_FLUXNET_proc(self):
                
        ec_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                     transforms = hourly_to_daily,
                                     NEE_partitioning_method = None,
                                     version    = self.version, 
                                     site       = self.site)
        ec_data = []
        for var_ in list(self.vars.keys()):
            ec_vars = []
            ec_vars.append(Variable(self.vars[var_]['sourceVariableName'], depth = "shallow"))   
            ec_vars.append(Variable(self.vars[var_]['sourceVariableName'], depth = "deep"))   
            ds_ = ec_prov.get_data(ec_vars)    
            ds_ = xr.DataArray(ds_.to_array().values.reshape(len(list(ds_.keys())),-1,1,1), 
                         dims=['depth_FLUXNET', 'time', 'lat', 'lon'],
                         coords={'depth_FLUXNET': np.arange(len(list(ds_.keys()))) +1,
                                 'time': ds_.time.values,
                                 'lat': [self.lat],
                                 'lon': [self.lon]}).to_dataset(name=var_)
            ec_data.append(ds_)
        ec_data = xr.merge(ec_data)        
               
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            ec_data[var_] = ec_data[var_] * unit_scalar
            #if (np.nanmin(ec_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(ec_data[var_]) > self.vars[var_]['bounds'][1]):
            #    print("Values are out of bounds:\n"+var_)
                        
        return ec_data#.sel(time = slice(self.start_date, self.end_date))
    
class Climate_FLUXNET:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.version = config["FLUXNET_version"]
        self.vars = config["variables"]["climate_vars"]['FLUXNET_climate']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        
    def Climate_FLUXNET_proc(self):
        
        #%% Process data - DayMean
        DayMeanClimate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = self.cubepath,
                                                                         transforms     = hourly_to_daily,
                                                                         version        = self.version, 
                                                                         site           =self.site)
        DayMeanClimate_FLUXNET_vars = []
        for var_ in [x for x in list(self.vars.keys()) if 'DayMean' in x]:
           DayMeanClimate_FLUXNET_vars.append(Variable(self.vars[var_]['sourceVariableName'], unit = self.vars[var_]['sourceVariableUnit']))   
        DayMeanClimate_FLUXNET_data = DayMeanClimate_FLUXNET_prov.get_data(DayMeanClimate_FLUXNET_vars)
        
        #%% Rename data
        DayMeanClimate_FLUXNET_data = DayMeanClimate_FLUXNET_data.rename(**{v:v+'_DayMean_FLUXNET' for v in list(set(DayMeanClimate_FLUXNET_data.variables)-set(DayMeanClimate_FLUXNET_data.coords))})
                
        #%% Process data - DayTime
        DayTimeClimate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                                                         version     = self.version, 
                                                                         site        = self.site)
        daytime_hourly_to_daily = [transformers.DaytimeMean(DayTimeClimate_FLUXNET_prov.get_data(Variable('SW_IN_POT')), 
                                                            DayTimeClimate_FLUXNET_prov.get_data(Variable('SW_IN')))]    
        DayTimeClimate_FLUXNET_prov.add_transform(daytime_hourly_to_daily)
    
        DayTimeClimate_FLUXNET_vars = []        
        for var_ in [x for x in list(self.vars.keys()) if 'DayTime' in x]:
           DayTimeClimate_FLUXNET_vars.append(Variable(self.vars[var_]['sourceVariableName'], unit = self.vars[var_]['sourceVariableUnit']))   
        DayTimeClimate_FLUXNET_data = DayTimeClimate_FLUXNET_prov.get_data(DayTimeClimate_FLUXNET_vars)
        
        #%% Rename data   
        DayTimeClimate_FLUXNET_data = DayTimeClimate_FLUXNET_data.rename(**{v:v+'_DayTime_FLUXNET' for v in list(set(DayTimeClimate_FLUXNET_data.variables)-set(DayTimeClimate_FLUXNET_data.coords))})
        
        #%% Merge DayMean and DayTime variable
        Climate_FLUXNET_data = xr.merge([DayMeanClimate_FLUXNET_data, DayTimeClimate_FLUXNET_data])
                
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            Climate_FLUXNET_data[var_] = Climate_FLUXNET_data[var_] * unit_scalar
            #if (np.nanmin(Climate_FLUXNET_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(Climate_FLUXNET_data[var_]) > self.vars[var_]['bounds'][1]):
            #    print("Values are out of bounds:\n"+var_)
                
        #%% Restructure data
        Climate_FLUXNET_data = data_structure_temporal_NoDepth(Climate_FLUXNET_data, self.lat, self.lon)        
        
        return Climate_FLUXNET_data#.sel(time = slice(self.start_date, self.end_date))
    
class Climate_ERA5:   
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.version = config["FLUXNET_version"]
        self.vars = config["variables"]["climate_vars"]["ERA5_climate"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
    
    def Climate_ERA5_proc(self):
        
        #%% Process data - DayMean
        DayMeanClimate_ERA5_prov    = era5.ERA5Provider(cubepath = self.cubepath,
                                                        transforms     = hourly_to_daily,
                                                        version        = self.version, 
                                                        site           =self.site)
        DayMeanClimate_ERA5_vars = []
        for var_ in [x for x in list(self.vars.keys()) if 'DayMean' in x]:
           DayMeanClimate_ERA5_vars.append(Variable(self.vars[var_]['sourceVariableName'], unit = self.vars[var_]['sourceVariableUnit']))   
        DayMeanClimate_ERA5_data = DayMeanClimate_ERA5_prov.get_data(DayMeanClimate_ERA5_vars)
        
        #%% Rename data
        DayMeanClimate_ERA5_data = DayMeanClimate_ERA5_data.rename(**{v:v+'_DayMean_ERA5' for v in list(set(DayMeanClimate_ERA5_data.variables)-set(DayMeanClimate_ERA5_data.coords))})
                
        #%% Process data - DayTime
        DayTimeClimate_ERA5_prov    = era5.ERA5Provider(cubepath   = self.cubepath,
                                                        version    = self.version, 
                                                        site       = self.site)
        daytime_hourly_to_daily = [transformers.DaytimeMean(DayTimeClimate_ERA5_prov.get_data(Variable('SW_IN_POT')), 
                                                            DayTimeClimate_ERA5_prov.get_data(Variable('SW_IN')))]    
        DayTimeClimate_ERA5_prov.add_transform(daytime_hourly_to_daily)
    
        DayTimeClimate_ERA5_vars = []        
        for var_ in [x for x in list(self.vars.keys()) if 'DayTime' in x]:
           DayTimeClimate_ERA5_vars.append(Variable(self.vars[var_]['sourceVariableName'], unit = self.vars[var_]['sourceVariableUnit']))   
        DayTimeClimate_ERA5_data = DayTimeClimate_ERA5_prov.get_data(DayTimeClimate_ERA5_vars)
        
        #%% Rename data   
        DayTimeClimate_ERA5_data = DayTimeClimate_ERA5_data.rename(**{v:v+'_DayTime_ERA5' for v in list(set(DayTimeClimate_ERA5_data.variables)-set(DayTimeClimate_ERA5_data.coords))})
        
        #%% Merge DayMean and DayTime variable
        Climate_ERA5_data = xr.merge([DayMeanClimate_ERA5_data, DayTimeClimate_ERA5_data])
                
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            Climate_ERA5_data[var_] = Climate_ERA5_data[var_] * unit_scalar
            if (np.nanmin(Climate_ERA5_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(Climate_ERA5_data[var_]) > self.vars[var_]['bounds'][1]):
                print("Values are out of bounds:\n"+var_)
        
        #%% Restructure data
        Climate_ERA5_data = data_structure_temporal_NoDepth(Climate_ERA5_data, self.lat, self.lon)        
                
        return Climate_ERA5_data#.sel(time = slice(self.start_date, self.end_date))

class gapfill_climate:
    
        def __init__(self,Climate_site_data, Climate_downscaled_data, config):    
            self.Climate_site_data = Climate_site_data.copy()
            self.Climate_downscaled_data = Climate_downscaled_data.copy()
            self.QC_climate = config["QC_climate"]            
            
        def gapfill_climate_proc(self):
            
            climate_FLUXNET_gapfilled = []
            for var_ in list(self.Climate_downscaled_data.keys()):                
                var_FLUXNET = var_.replace('ERA5','FLUXNET')
                try:
                    var_FLUXNET_afterQA = self.Climate_site_data[var_FLUXNET].where(self.Climate_site_data[var_FLUXNET.replace('Day','QC_Day')] > self.QC_climate)
                except:
                    var_FLUXNET_afterQA = self.Climate_site_data[var_FLUXNET]
                    print(f"No QC variable for {var_FLUXNET} - No QC applied")
                
                #%% Perform the gapfilling
                climate_FLUXNET_gapfilled.append(var_FLUXNET_afterQA.where(np.isfinite(var_FLUXNET_afterQA), self.Climate_downscaled_data[var_]).to_dataset(name = var_FLUXNET))
                
            #%% Rename variable
            climate_FLUXNET_gapfilled = xr.merge(climate_FLUXNET_gapfilled)
            climate_FLUXNET_gapfilled = climate_FLUXNET_gapfilled.rename(**{v:v+'_gapfilled' for v in list(set(climate_FLUXNET_gapfilled.variables)-set(climate_FLUXNET_gapfilled.coords))})   
            
            return climate_FLUXNET_gapfilled
    
class MODIS_MCD43A:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.vars = config["variables"]["MODIS_vars"]["MODIS_MCD43A"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
    
    def MODIS_MCD43A_proc(self):
        
        #%% Process data 
        modis_prov = modis_MCD43A.MCD43A(cubepath =  self.cubepath,  site =self.site)
        modis_vars = []
        for var_ in list(self.vars.keys()):
            modis_vars.append(Variable(self.vars[var_]['sourceVariableName'], unit = self.vars[var_]['sourceVariableUnit']))   
        modis_data = modis_prov.get_data(modis_vars)
        
        #%% Rename data
        modis_data = modis_data.rename(**{v:v+'_MCD43A' for v in list(set(modis_data.variables)-set(modis_data.coords))})
        
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            modis_data[var_] = modis_data[var_] * unit_scalar
            if (np.nanmin(modis_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(modis_data[var_]) > self.vars[var_]['bounds'][1]):
                print("Values are out of bounds:\n"+var_)
        
        #%% Restructure data
        modis_data = data_structure_temporal_NoDepth(modis_data, self.lat, self.lon)               
                
        return modis_data#.sel(time = slice(self.start_date, self.end_date))
    
class MODIS_MxD11A:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.vars = config["variables"]["MODIS_vars"]["MODIS_MxD11A"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
    
    def MODIS_MxD11A_proc(self):
        
        #%% Process data 
        modis_prov = modis_MxD11A.MxD11A(cubepath =  self.cubepath, 
                                         satellite=None, day_night=None,
                                         site =self.site)
        LST_Day_AQUA = modis_prov.get_data(Variable('LST', satellite='AQUA', day_night='Day'))
        LST_Night_AQUA = modis_prov.get_data(Variable('LST', satellite='AQUA', day_night='Night'))
        LST_Day_TERRA = modis_prov.get_data(Variable('LST', satellite='TERRA', day_night='Day'))
        LST_Night_TERRA = modis_prov.get_data(Variable('LST', satellite='TERRA', day_night='Night'))
        
        for var_ in list(self.vars.keys()):
        
            if var_ == 'LST_DayMean_AQUA_MxD11A':
                LST_DayMean_AQUA = xr.merge([LST_Day_AQUA, LST_Night_AQUA]).to_array().mean(axis=0).to_dataset(name = 'LST_DayMean_AQUA_MxD11A')
            
            elif var_ == 'LST_DayMean_TERRA_MxD11A':
                LST_DayMean_TERRA = xr.merge([LST_Day_TERRA, LST_Night_TERRA]).to_array().mean(axis=0).to_dataset(name = 'LST_DayMean_TERRA_MxD11A')
            
            elif var_ == 'LST_DayTime_AQUA_MxD11A':
                LST_DayTime_AQUA = LST_Day_AQUA.to_dataset(name = 'LST_DayTime_AQUA_MxD11A')
            
            elif var_ == 'LST_DayTime_TERRA_MxD11A':
                LST_DayTime_TERRA = LST_Day_TERRA.to_dataset(name = 'LST_DayTime_TERRA_MxD11A')                    
                    
        modis_data = xr.merge([LST_DayMean_AQUA, LST_DayMean_TERRA, LST_DayTime_AQUA, LST_DayTime_TERRA])
        
        #%% Restructure data
        modis_data = data_structure_temporal_NoDepth(modis_data, self.lat, self.lon)   
        
        return modis_data#.sel(time = slice(self.start_date, self.end_date))

class EvapTrans_Nelson2018:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]        
        self.datapath = config["DataPath"]["TEA_data"] + '/' + self.site + '/TEA/L1_daily/'
        self.vars = config["variables"]["EvapTrans_vars"]["Nelson2018_EvapTrans"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
            
    def EvapTrans_Nelson2018_proc(self):
        
        #%% Process data
        TEA_data = []
        for var_ in list(self.vars.keys()):
            try:
                TEA_data.append(xr.open_mfdataset(self.datapath + "*.nc")[self.vars[var_]['sourceVariableName']].to_dataset().rename({self.vars[var_]['sourceVariableName']: var_}))
            except:
                TEA_data.append(xr.DataArray(np.nan, dims=['time'],coords={'time': [np.datetime64(str(self.start_date))]}).to_dataset(name = var_))
        TEA_data = xr.merge(TEA_data)
        
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            TEA_data[var_] = TEA_data[var_] * unit_scalar
            if (np.nanmin(TEA_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(TEA_data[var_]) > self.vars[var_]['bounds'][1]):
                print("Values are out of bounds:\n"+var_)
        
        #%% Restructure data
        TEA_data = data_structure_temporal_NoDepth(TEA_data, self.lat, self.lon)        
                        
        return TEA_data#.sel(time = slice(self.start_date, self.end_date))
    
class insitu_AGB:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]                
        self.datapath = config["DataPath"]["insituAGB"]
        self.vars = config["variables"]['biomass_vars']["AG_BIOMASS_TT_vars"]['AG_BIOMASS_TT_insitu']['sourceVariableName']
        self.scalar =  config["variables"]['biomass_vars']["AG_BIOMASS_TT_vars"]['AG_BIOMASS_TT_insitu']['source2sindbadUnit']
        self.bounds = config["variables"]['biomass_vars']["AG_BIOMASS_TT_vars"]['AG_BIOMASS_TT_insitu']['bounds']
    
    def insitu_AGB_proc(self):
        
        #%% Process data 
        insitu_AGB_df = pd.read_csv(self.datapath)
        insitu_AGB_data = insitu_AGB_df.loc[insitu_AGB_df['Site_ID'] == self.site]
        if insitu_AGB_data.size == 0:
            insitu_AGB_data = xr.DataArray(np.nan, dims=['time'],
                                            coords={'time': [np.datetime64(str(self.start_date))]}).to_dataset(name = self.vars)
        else:
            insitu_AGB_data = xr.DataArray(insitu_AGB_data['insitu_AGB'].values, dims=['time'],
                                           coords={'time':  [np.datetime64(str(y)) for y in insitu_AGB_data.year.values]}).to_dataset(name = self.vars)
            
        #%% Apply unit correction and check value bounds
        unit_scalar = int(self.scalar)
        insitu_AGB_data[self.vars] = insitu_AGB_data[self.vars] * unit_scalar
        if (np.nanmin(insitu_AGB_data[self.vars]) < self.bounds[0]) | (np.nanmax(insitu_AGB_data[self.vars]) > self.bounds[1]):
            print("Values are out of bounds:\n"+self.vars)
            
        #%% Restructure data
        insitu_AGB_data = data_structure_temporal_NoDepth(insitu_AGB_data, self.lat, self.lon)        
                        
        return insitu_AGB_data
    
class gloBbiomass_AGB:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']     
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]           
        self.datapath = config["DataPath"]["GlobBiomassAGB"]
        self.vars = config["variables"]['biomass_vars']["AG_BIOMASS_TT_vars"]['AG_BIOMASS_TT_gloBbiomass']['sourceVariableName']
        self.scalar =  config["variables"]['biomass_vars']["AG_BIOMASS_TT_vars"]['AG_BIOMASS_TT_gloBbiomass']['source2sindbadUnit']
        self.bounds = config["variables"]['biomass_vars']["AG_BIOMASS_TT_vars"]['AG_BIOMASS_TT_gloBbiomass']['bounds']
    
    def gloBbiomass_AGB_proc(self):
        
        #%% Process data     
        gloBbiomass_AGB_df = pd.read_csv(self.datapath)
        gloBbiomass_AGB_data = np.array(gloBbiomass_AGB_df.loc[gloBbiomass_AGB_df['Site_ID'] == self.site].drop(['Site_ID', 'lon', 'lat'], axis=1), dtype=float)
        if gloBbiomass_AGB_data.size == 0:
            gloBbiomass_AGB_data = xr.DataArray(np.nan, dims=['time'],
                                       coords={'time': [np.datetime64(str(self.start_date))]}).to_dataset(name = self.vars)
        else:
            gloBbiomass_AGB_data[gloBbiomass_AGB_data == -9999] = np.nan
            gloBbiomass_AGB_data = np.nanmedian(gloBbiomass_AGB_data)
            gloBbiomass_AGB_data = xr.DataArray(gloBbiomass_AGB_data, dims=['time'],
                                       coords={'time': [np.datetime64(str(2010))]}).to_dataset(name = self.vars)
        
        #%% Apply unit correction and check value bounds
        unit_scalar = int(self.scalar)
        gloBbiomass_AGB_data[self.vars] = gloBbiomass_AGB_data[self.vars] * unit_scalar
        if (np.nanmin(gloBbiomass_AGB_data[self.vars]) < self.bounds[0]) | (np.nanmax(gloBbiomass_AGB_data[self.vars]) > self.bounds[1]):
            print("Values are out of bounds:\n"+self.vars)        
        
        #%% Restructure data
        gloBbiomass_AGB_data = data_structure_temporal_NoDepth(gloBbiomass_AGB_data, self.lat, self.lon)        
                
        return gloBbiomass_AGB_data
    
class hansen_treecover:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]                
        self.datapath = config["DataPath"]["HansenTreecover"]
        self.sourceVariableName = config["variables"]["TREECOVER_vars"]['TREECOVER_RS']['TREECOVER_hansen']['sourceVariableName']
        self.vars = config["variables"]["TREECOVER_vars"]['TREECOVER_RS']['TREECOVER_hansen']['nameShort']        
        self.scalar =  config["variables"]["TREECOVER_vars"]['TREECOVER_RS']['TREECOVER_hansen']['source2sindbadUnit']
        self.bounds = config["variables"]["TREECOVER_vars"]['TREECOVER_RS']['TREECOVER_hansen']['bounds']
    
    def hansen_treecover_proc(self):
        
        #%% Process data     
        hansen_treecover_df = pd.read_csv(self.datapath)
        hansen_treecover_data = np.array(hansen_treecover_df.loc[hansen_treecover_df['SiteID'] == self.site][self.sourceVariableName], dtype=float)
        if hansen_treecover_data.size == 0:
            hansen_treecover_data = xr.DataArray(np.nan, dims=['time'],
                                       coords={'time': [np.datetime64(str(self.start_date))]}).to_dataset(name = self.vars)
        else:
            hansen_treecover_data = xr.DataArray(hansen_treecover_data, dims=['time'],
                                       coords={'time': [np.datetime64(str(2010))]}).to_dataset(name = self.vars)
        
        #%% Apply unit correction and check value bounds
        unit_scalar = int(self.scalar)
        hansen_treecover_data[self.vars] = hansen_treecover_data[self.vars] * unit_scalar
        if (np.nanmin(hansen_treecover_data[self.vars]) < self.bounds[0]) | (np.nanmax(hansen_treecover_data[self.vars]) > self.bounds[1]):
            print("Values are out of bounds:\n"+self.vars)        
        
        #%% Restructure data
        hansen_treecover_data = data_structure_temporal_NoDepth(hansen_treecover_data, self.lat, self.lon)        
        
        return hansen_treecover_data

class forest_age:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.datapath = config["DataPath"]["forest_age"]
        self.vars = config["variables"]['Forest_age_vars']["DIST_FOREST_FRAC"]['DIST_FOREST_FRAC_Besnard2018']['sourceVariableName']
        self.scalar =  config["variables"]['Forest_age_vars']["DIST_FOREST_FRAC"]['DIST_FOREST_FRAC_Besnard2018']['source2sindbadUnit']
        self.bounds = config["variables"]['Forest_age_vars']["DIST_FOREST_FRAC"]['DIST_FOREST_FRAC_Besnard2018']['bounds']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        
    def forest_age_proc(self):
        
        #%% Process data 
        forest_age_df = pd.read_csv(self.datapath)
        
        LAST_DISTURBANCE_DATE =  forest_age_df.loc[forest_age_df['Site_ID'] == self.site]['Plantation_Date_max'].values.astype(np.datetime64)
        if LAST_DISTURBANCE_DATE.size ==0:
            LAST_DISTURBANCE_DATE = [np.nan]
            forest_age_data = xr.DataArray(np.nan, dims=['time'], coords={'time': [np.datetime64(str(self.start_date))]}).to_dataset(name = self.vars)
        else:    
            date_ = np.arange(np.datetime64(self.start_date), np.datetime64(self.end_date) + np.timedelta64(1,'D') ,dtype="M8[D]")        
            forest_age_data = np.concatenate([(day_ - LAST_DISTURBANCE_DATE).astype('int') / 365.25 for day_ in date_] )
            forest_age_data = xr.DataArray(forest_age_data, dims=['time'],
                                       coords={'time': date_})
            forest_age_data = forest_age_data.where(forest_age_data==0)
            forest_age_data = forest_age_data.where(np.isfinite(forest_age_data),1)
            forest_age_data = forest_age_data.to_dataset(name = self.vars)
           
        #%% Apply unit correction and check value bounds
        unit_scalar = int(self.scalar)
        forest_age_data[self.vars] = forest_age_data[self.vars] * unit_scalar
        if (np.nanmin(forest_age_data[self.vars]) < self.bounds[0]) | (np.nanmax(forest_age_data[self.vars]) > self.bounds[1]):
            print("Values are out of bounds:\n"+self.vars)
            
        #%% Restructure data
        forest_age_data = data_structure_temporal_NoDepth(forest_age_data, self.lat, self.lon)        
        
        return forest_age_data, LAST_DISTURBANCE_DATE
    
class insitu_VegFrac:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.datapath = config["DataPath"]["PFT"]
        self.vars = config["variables"]["Veg_Frac_vars"]['Veg_Frac_insitu']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]        
        
    def insitu_VegFrac_proc(self):
        
        #%% Process data     
        PFT_df = pd.read_csv(self.datapath)
        PFT_ = PFT_df.loc[PFT_df['Site_ID'] == self.site]
        date_ = np.arange(np.datetime64(self.start_date), np.datetime64(self.end_date) + np.timedelta64(1,'D') ,dtype="M8[D]")        
        if PFT_.size == 0:
            PFT_ = np.nan
            Frac_tree_insitu = xr.Dataset({self.vars['Frac_tree_insitu']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_grass_insitu = xr.Dataset({self.vars['Frac_grass_insitu']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_shrub_insitu = xr.Dataset({self.vars['Frac_shrub_insitu']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_crop_insitu = xr.Dataset({self.vars['Frac_crop_insitu']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_savanna_insitu = xr.Dataset({self.vars['Frac_savanna_insitu']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_wetland_insitu = xr.Dataset({self.vars['Frac_wetland_insitu']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
        
        else:
            PFT_ = PFT_['PFT'].values
            
            #%% Tree fraction
            if PFT_[0] in 'DBF|MF|ENF|EBF':
                Frac_tree_insitu = xr.Dataset({self.vars['Frac_tree_insitu']['sourceVariableName']: xr.DataArray(data = 1, dims = ['time'], coords = {'time': date_})})
            else:
                Frac_tree_insitu = xr.Dataset({self.vars['Frac_tree_insitu']['sourceVariableName']: xr.DataArray(data =0, dims = ['time'], coords = {'time': date_})})
            
            #%% Grass fraction
            if PFT_[0] in 'GRA':
                Frac_grass_insitu = xr.Dataset({self.vars['Frac_grass_insitu']['sourceVariableName']: xr.DataArray(data = 1, dims = ['time'], coords = {'time': date_})})
            else:
                Frac_grass_insitu = xr.Dataset({self.vars['Frac_grass_insitu']['sourceVariableName']: xr.DataArray(data = 0, dims = ['time'], coords = {'time': date_})})
            
            #%% shrub fraction
            if PFT_[0] in 'OSH':
                Frac_shrub_insitu = xr.Dataset({self.vars['Frac_shrub_insitu']['sourceVariableName']: xr.DataArray(data = 1, dims = ['time'], coords = {'time': date_})})
            else:
                Frac_shrub_insitu = xr.Dataset({self.vars['Frac_shrub_insitu']['sourceVariableName']: xr.DataArray(data = 0, dims = ['time'], coords = {'time': date_})})
            
            #%% Crop fraction
            if PFT_[0] in 'CRO':
                Frac_crop_insitu = xr.Dataset({self.vars['Frac_crop_insitu']['sourceVariableName']: xr.DataArray(data = 1, dims = ['time'], coords = {'time': date_})})
            else:
                Frac_crop_insitu = xr.Dataset({self.vars['Frac_crop_insitu']['sourceVariableName']: xr.DataArray(data = 0, dims = ['time'], coords = {'time': date_})})
            
            #%% Savanna fraction
            if PFT_[0] in 'SAV|WSA':
                Frac_savanna_insitu = xr.Dataset({self.vars['Frac_savanna_insitu']['sourceVariableName']: xr.DataArray(data = 1, dims = ['time'], coords = {'time': date_})})
            else:
                Frac_savanna_insitu = xr.Dataset({self.vars['Frac_savanna_insitu']['sourceVariableName']: xr.DataArray(data = 0, dims = ['time'], coords = {'time': date_})})
            
            #%% Wetland fraction
            if PFT_[0] in 'WET':
                Frac_wetland_insitu = xr.Dataset({self.vars['Frac_wetland_insitu']['sourceVariableName']: xr.DataArray(data = 1, dims = ['time'], coords = {'time': date_})})
            else:
                Frac_wetland_insitu = xr.Dataset({self.vars['Frac_wetland_insitu']['sourceVariableName']: xr.DataArray(data = 0, dims = ['time'], coords = {'time': date_})})
            
        #%% Merge class
        Frac_Veg_insitu = xr.merge([Frac_tree_insitu, Frac_grass_insitu, Frac_shrub_insitu, Frac_crop_insitu, Frac_savanna_insitu, Frac_wetland_insitu])
        
        #%% Restructure data
        Frac_Veg_insitu = data_structure_temporal_NoDepth(Frac_Veg_insitu, self.lat, self.lon)   
        
        return PFT_, Frac_Veg_insitu, 

class hilda_VegFrac:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.datapath = config["DataPath"]["hilda"]
        self.vars = config["variables"]["Veg_Frac_vars"]['Veg_Frac_hilda']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]        
        
    def hilda_VegFrac_proc(self):
        
        #%% Process data     
        hilda_data = xr.open_dataset(self.datapath).sel(time = slice(int(self.start_date[0:4]),int(self.end_date[0:4])))
        hilda_site = hilda_data.sel(latitude = self.lat, longitude = self.lon, method='nearest').LULC_states
        hilda_site = np.round(hilda_site, -1)        
        date_ = np.arange(np.datetime64(self.start_date), np.datetime64(self.end_date) + np.timedelta64(1,'D') ,dtype="M8[D]")        
        if hilda_site.size == 0:
            Frac_tree_hilda = xr.Dataset({self.vars['Frac_tree_hilda']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_grass_hilda = xr.Dataset({self.vars['Frac_grass_hilda']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_shrub_hilda = xr.Dataset({self.vars['Frac_shrub_hilda']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_crop_hilda = xr.Dataset({self.vars['Frac_crop_hilda']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            Frac_other_hilda = xr.Dataset({self.vars['Frac_other_hilda']['sourceVariableName']: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
        
        else:
            Frac_crop_hilda = hilda_site.where(hilda_site == 20, 0)
            Frac_crop_hilda[Frac_crop_hilda>0] = 1
            Frac_tree_hilda = hilda_site.where(hilda_site == 40, 0)
            Frac_tree_hilda[Frac_tree_hilda>0] = 1            
            Frac_grass_hilda = hilda_site.where(hilda_site == 30, 0)
            Frac_grass_hilda[Frac_grass_hilda>0] = 1
            Frac_shrub_hilda = hilda_site.where(hilda_site == 60, 0)
            Frac_shrub_hilda[Frac_shrub_hilda>0] = 1
            Frac_other_hilda = hilda_site.where(hilda_site == 70, 0)
            Frac_other_hilda[Frac_other_hilda>0] = 1
            annual_crop = []
            annual_tree = []
            annual_grass = []
            annual_shrub = []
            annual_other = []            
            for year_ in hilda_site.time:
                date_ = np.arange(str(int(year_.values)) + '-01' + '-01', str(int(year_.values) +1) + '-01' + '-01'  ,dtype="M8[D]")
                annual_crop.append(xr.Dataset({self.vars['Frac_crop_hilda']['sourceVariableName']: xr.DataArray(data = int(Frac_crop_hilda.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_})}))
                annual_tree.append(xr.Dataset({self.vars['Frac_tree_hilda']['sourceVariableName']: xr.DataArray(data = int(Frac_tree_hilda.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_})}))
                annual_grass.append(xr.Dataset({self.vars['Frac_grass_hilda']['sourceVariableName']: xr.DataArray(data = int(Frac_grass_hilda.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_})}))
                annual_shrub.append(xr.Dataset({self.vars['Frac_shrub_hilda']['sourceVariableName']: xr.DataArray(data = int(Frac_shrub_hilda.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_})}))
                annual_other.append(xr.Dataset({self.vars['Frac_other_hilda']['sourceVariableName']: xr.DataArray(data = int(Frac_other_hilda.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_})}))
            Frac_crop_hilda = xr.merge(annual_crop)   
            Frac_tree_hilda = xr.merge(annual_tree)   
            Frac_grass_hilda = xr.merge(annual_grass)   
            Frac_shrub_hilda = xr.merge(annual_shrub)   
            Frac_other_hilda = xr.merge(annual_other)

        #%% Merge class
        Frac_Veg_hilda = xr.merge([Frac_tree_hilda, Frac_grass_hilda, Frac_shrub_hilda, Frac_crop_hilda, Frac_other_hilda])
        
        #%% Restructure data
        Frac_Veg_hilda = data_structure_temporal_NoDepth(Frac_Veg_hilda, self.lat, self.lon)
             
        return Frac_Veg_hilda
    
class soilGrids_soilTexture:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.datapath = config["DataPath"]["SoilGrids"]
        self.vars = config["variables"]["Soil_texture_vars"]['Soil_texture_SoilGrids']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
                
    def soilGrids_soilTexture_proc(self):
        
        #%% Process data     
        soilGrids_soilTexture_df = pd.read_csv(self.datapath)
        soilGrids_vars = []
        for var_ in list(self.vars.keys()):
            soilTexture_var = soilGrids_soilTexture_df.filter(like = self.vars[var_]['sourceVariableName'], axis=1).filter(regex='^((?!Sync).)*$').loc[soilGrids_soilTexture_df['SiteID'] == self.site]
            soilGrids_soilTexture_data = soilTexture_var.reindex(sorted(soilTexture_var.columns), axis=1).values.reshape(-1).astype('d')
            
            if soilGrids_soilTexture_data.size == 0:
                soilGrids_soilTexture_data = xr.DataArray(np.repeat(np.nan, len(soilTexture_var.columns)).reshape(len(soilTexture_var.columns), 1, 1), 
                                                          dims=['depth_soilGrids','lat', 'lon'],
                                                          coords={'depth_soilGrids': np.arange(len(soilTexture_var.columns))+1,
                                                                  'lat': [self.lat],
                                                                  'lon': [self.lon]}).to_dataset(name = var_)
            else:
                soilGrids_soilTexture_data = xr.DataArray(soilGrids_soilTexture_data.reshape(soilGrids_soilTexture_data.size, 1, 1), 
                                                          dims=['depth_soilGrids','lat', 'lon'],
                                                          coords={'depth_soilGrids': np.arange(soilGrids_soilTexture_data.size)+1,
                                                                  'lat': [self.lat],
                                                                  'lon': [self.lon]}).to_dataset(name = var_)
            
            #%% Apply unit correction and check value bounds
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            soilGrids_soilTexture_data[self.vars[var_]['nameShort']] = soilGrids_soilTexture_data[self.vars[var_]['nameShort']] * unit_scalar
            if (np.nanmin(soilGrids_soilTexture_data[self.vars[var_]['nameShort']]) < self.vars[var_]['bounds'][0]) | (np.nanmax(soilGrids_soilTexture_data[self.vars[var_]['nameShort']]) > self.vars[var_]['bounds'][1]):
                #soilGrids_soilTexture_data[var_] = np.repeat(np.nan, soilGrids_soilTexture_data.depth_soilGrids.shape[0]).reshape(soilGrids_soilTexture_data.depth_soilGrids.shape[0], 1, 1)
                print(f'Values are out of bounds for {var_}')
            soilGrids_vars.append(soilGrids_soilTexture_data)
        
        soilGrids_vars = xr.merge(soilGrids_vars)        
        return soilGrids_vars
    
class data_harm:
    def __init__(self, data, config):
        self.data = data
        self.config = config
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        
    def data_harm_proc(self):
        
        #%Sort variable alphabeticatly
        var_ = sorted(list(self.data.keys()))
        self.data = self.data[var_]    
        
        #%% Add metadata
        for dat_type in self.config["variables"].keys():
                               
                for source_ in list(self.config["variables"][dat_type].keys()):
                    
                    for var_ in list(self.config["variables"][dat_type][source_].keys()):
                        if 'ERA5' in var_:
                            self.data[var_.replace('ERA5',"FLUXNET_gapfilled")] = self.data[var_.replace('ERA5',"FLUXNET_gapfilled")].assign_attrs(bounds = str([np.nanmin(self.data[var_.replace('ERA5',"FLUXNET_gapfilled")]), np.nanmax(self.data[var_.replace('ERA5',"FLUXNET_gapfilled")])]),
                                                                                                                                       QA_applied = self.config['QC_climate'],
                                                                                                                                       units = self.config["variables"][dat_type][source_][var_]['variableUnit'],
                                                                                                                                       short_name = var_.replace('ERA5',"FLUXNET_gapfilled"),
                                                                                                                                       long_name = var_.replace('ERA5',"FLUXNET_gapfilled") + ' - FLUXNET climate data gapfilled with ERA5',
                                                                                                                                       sourceDataProductName = self.config["variables"][dat_type][source_][var_]['sourceDataProductName'],
                                                                                                                                       publication = self.config["variables"][dat_type][source_][var_]['publication'],
                                                                                                                                       dataPath = self.config["variables"][dat_type][source_][var_]['dataPath'],
                                                                                                                                       )
                            self.data[var_] = self.data[var_].assign_attrs(bounds = str([np.nanmin(self.data[var_]), np.nanmax(self.data[var_])]),
                                                               units = self.config["variables"][dat_type][source_][var_]['variableUnit'],
                                                               short_name = self.config["variables"][dat_type][source_][var_]['nameShort'],
                                                               long_name = self.config["variables"][dat_type][source_][var_]['nameLong'],
                                                               sourceDataProductName = self.config["variables"][dat_type][source_][var_]['sourceDataProductName'],
                                                               publication = self.config["variables"][dat_type][source_][var_]['publication'],
                                                               dataPath = self.config["variables"][dat_type][source_][var_]['dataPath'])
                        else:
                            self.data[var_] = self.data[var_].assign_attrs(bounds = str([np.nanmin(self.data[var_]), np.nanmax(self.data[var_])]),
                                                               units = self.config["variables"][dat_type][source_][var_]['variableUnit'],
                                                               short_name = self.config["variables"][dat_type][source_][var_]['nameShort'],
                                                               long_name = self.config["variables"][dat_type][source_][var_]['nameLong'],
                                                               sourceDataProductName = self.config["variables"][dat_type][source_][var_]['sourceDataProductName'],
                                                               publication = self.config["variables"][dat_type][source_][var_]['publication'],
                                                               dataPath = self.config["variables"][dat_type][source_][var_]['dataPath'])
                     
        return self.data.sel(time = slice(self.start_date, self.end_date))
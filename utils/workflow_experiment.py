#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
#%% Load library
from numpy.lib.ufunclike import _deprecate_out_named_y
import xarray as xr
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, date
import importlib
import logging
sys.path.append(os.path.join(os.getcwd(), '../'))
import copy
from extractors.cliff_gapfill import cliff_gapfill
import collections
import utils.shared_utils as shut


def get_site_info(site):
    import fluxcom as flx
    lat = flx.site_cube.site_lat_lon.lat[site]
    lon = flx.site_cube.site_lat_lon.lon[site]
    site_info = {"site_ID": site, "latitude": lat, "longitude": lon}
    return site_info


def load_extractor(mod_name):
    mod = importlib.import_module("extractors." + mod_name)
    return mod


def add_date_variables(dat):
    julian_base_date = date(1592, 10, 15)
    julian_base_date_str = julian_base_date.strftime("%Y-%m-%d")
    time = dat.time
    dates = pd.DatetimeIndex(time.values)
    year = dates.year.to_numpy()
    month = dates.month.to_numpy()
    day = dates.day.to_numpy()
    hour = dates.hour.to_numpy()  #+ 1

    julian_day = np.zeros_like(year)
    for _ind in range(len(year)):
        data_date = date(year[_ind], month[_ind], day[_ind])
        delta = data_date - julian_base_date
        julian_day[_ind] = delta.days

    x_year = xr.Variable(["time"],
                         year,
                         attrs={
                             'long_name': 'year',
                             '_FillValue': -9999
                         })
    dat = dat.assign(year=x_year)
    x_month = xr.Variable(["time"],
                          month,
                          attrs={
                              'long_name': 'month of year',
                              '_FillValue': -9999
                          })
    dat = dat.assign(month=x_month)
    x_day = xr.Variable(["time"],
                        day,
                        attrs={
                            'long_name': 'day of year',
                            '_FillValue': -9999
                        })
    dat = dat.assign(day=x_day)
    x_hour = xr.Variable(["time"],
                         hour,
                         attrs={
                             'long_name': 'hour of day',
                             '_FillValue': -9999
                         })
    dat = dat.assign(hour=x_hour)
    x_julian_day = xr.Variable(["time"],
                               julian_day,
                               attrs={
                                   'long_name':
                                   f'days since {julian_base_date_str}',
                                   'FillValue': np.nan
                               })
    dat = dat.assign(julian_day=x_julian_day)
    return dat


def get_bounds_str(data, bounds_in):
    bounds_str = f'bounds: data: [{round(np.nanmin(data.values),3)}, {round(np.nanmax(data.values),3)}], json: [{bounds_in[0]}, {bounds_in[1]}]'
    return bounds_str


def plot_dianostic_figures(ds, site_info, exp_settings, site, resampled=None):
    import matplotlib.pyplot as plt
    if resampled is None:
        data_freq = exp_settings["temporal_resolution"]
    else:
        data_freq = resampled
    site_fig_dir = os.path.join(exp_settings['OutPath']['figs'], data_freq)

    # site_fig_dir = os.path.join(exp_settings['OutPath']['figs'], site)
    os.makedirs(site_fig_dir, exist_ok=True)
    skipvars = 'month day year hour julian_day'.split()
    plotvars = list(set(list(ds.keys())) - set(skipvars))
    for data_key in exp_settings["sel_datasets"]:
        data_set = exp_settings["dataset"][data_key]
        for var_, var_info in data_set['variables'].items():
            if var_ in plotvars:
                dat_ = ds[var_]
                shut.log_and_print('plot',
                                   site,
                                   dat_.sourceDataProductName,
                                   var_,
                                   dat_.units,
                                   data_freq,
                                   tar_label='variable',
                                   src_label='units')
                # print(f'{dat_.sourceDataProductName} : {data_freq} | {site} : {var_}')
                if ('depth_FLUXNET' in dat_.dims) & ('time' in dat_.dims):
                    fig, ax = plt.subplots(1,
                                           1,
                                           figsize=(8, 4),
                                           gridspec_kw={
                                               'wspace': 0.35,
                                               'hspace': 0.5
                                           })
                    for layer_ in dat_.depth_FLUXNET.values:
                        ax.plot(
                            dat_.sel(
                                depth_FLUXNET=layer_).time.values.reshape(-1),
                            dat_.sel(depth_FLUXNET=layer_).values.reshape(-1),
                            lw=0.6,
                            label='Layer ' + str(layer_))
                    plt.legend()
                    ax.set_xlabel('time', size=8)
                    ax.set_ylabel(data_freq + " " + var_ + ' [' + dat_.units +
                                  ']',
                                  size=8)
                elif ('time' in dat_.dims) & (
                        'depth_FLUXNET' not in dat_.dims) & ('depth_soilGrids'
                                                             not in dat_.dims):
                    fig, ax = plt.subplots(1,
                                           1,
                                           figsize=(8, 4),
                                           gridspec_kw={
                                               'wspace': 0.35,
                                               'hspace': 0.5
                                           })
                    if not var_.endswith("_gfld"):
                        ax.plot(dat_.time.values.reshape(-1),
                                dat_.values.reshape(-1),
                                lw=0.6)
                    else:
                        # print(dat_, dat_.sourceDataProductName)
                        ax.plot(dat_.time.values.reshape(-1),
                                dat_.values.reshape(-1),
                                lw=0.6,
                                label=dat_.sourceDataProductName)
                        dat_ori_ = ds[var_info['with_gaps']]
                        # print(dat_ori_)
                        ax.plot(dat_ori_.time.values.reshape(-1),
                                dat_ori_.values.reshape(-1),
                                lw=0.6,
                                label=dat_ori_.sourceDataProductName)
                        plt.legend()

                    ax.set_xlabel('time', size=8)
                    ax.set_ylabel(data_freq + " " + var_ + ' [' + dat_.units +
                                  ']',
                                  size=8)
                elif 'depth_soilGrids' in dat_.dims:
                    if 'time' in dat_.dims:
                        ## time axis is automatically generated in the monthly data are created using resampling
                        dat_ = dat_.mean(dim='time', keep_attrs=True)
                    fig, ax = plt.subplots(1,
                                           1,
                                           figsize=(8, 4),
                                           gridspec_kw={
                                               'wspace': 0.35,
                                               'hspace': 0.5
                                           })
                    ax.scatter(dat_.depth_soilGrids.values, dat_.values, s=10)
                    ax.set_xlabel('depth_soilGrids', size=8)
                    ax.set_ylabel(data_freq + " " + var_ + ' [' + dat_.units +
                                  ']',
                                  size=8)
                ax.tick_params(labelsize=6)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.title(
                    f'{site_info["site_ID"]}:: {data_key} in json: {dat_.sourceDataProductName} with {exp_settings["FLUXNET_version"]}\n {get_bounds_str(dat_,var_info["bounds"])}',
                    multialignment='center',
                    va='top',
                    size=8,
                    y=1.03)
                plt.savefig(
                    site_fig_dir + "/" + site_info["site_ID"] + '_' + var_ +
                    "." + exp_settings['start_date'].split('-')[0] + "-" +
                    str(int(exp_settings['end_date'].split('-')[0])) + '.' +
                    data_freq + '.png',
                    dpi=300)
            else:
                logging.warning(
                    f'Variable {var_} from {data_key} dataset was not found in processed data.'
                )
    shut.log_and_print_sep()
    return

def add_gapfill_settings(exp_settings,
                         fill_src=None,
                         fill_tar=None,
                         dataset_key=None,
                         qc_thres=None):
    exp_settings['sel_datasets'] = np.append(exp_settings['sel_datasets'],
                                             dataset_key)
    exp_settings['dataset'][dataset_key] = copy.deepcopy(
        exp_settings['dataset'][fill_tar])
    exp_settings['dataset'][dataset_key]['gfld_src'] = fill_src
    exp_settings['dataset'][dataset_key]['gfld_tar'] = fill_tar

    src_suffix = exp_settings['dataset'][fill_src]['var_suffix']
    if len(src_suffix.strip()) > 0:
        src_suffix = '_' + src_suffix

    tar_suffix = exp_settings['dataset'][fill_tar]['var_suffix']
    if len(tar_suffix.strip()) > 0:
        tar_suffix = '_' + tar_suffix
    tar_info = copy.deepcopy(exp_settings['dataset'][fill_tar]['variables'])
    src_info = copy.deepcopy(exp_settings['dataset'][fill_src]['variables'])
    _v_tars = list(exp_settings['dataset'][fill_tar]['variables'].keys())

    for _v_tar in _v_tars:
        _v_tar_gfld = _v_tar + src_suffix + '_gfld'
        _v_src = _v_tar_gfld.replace(tar_suffix, '').replace('_gfld', '')
        gfl_var_info = copy.deepcopy(tar_info[_v_tar])
        if _v_src in src_info.keys():
            gfl_var_info[
                'sourceDataProductName'] = f"{tar_info[_v_tar]['sourceDataProductName']} gapfld-w {src_info[_v_src]['sourceDataProductName']}"
            gfl_var_info[
                'publication'] = f"{tar_info[_v_tar]['publication']}; {src_info[_v_src]['publication']}"
            gfl_var_info[
                'description'] = f"{tar_info[_v_tar]['sourceDataProductName']}:: {tar_info[_v_tar]['description']};{src_info[_v_src]['sourceDataProductName']}:: {src_info[_v_src]['description']}"
            gfl_var_info[
                'sourceVariableName'] = f"{tar_info[_v_tar]['sourceDataProductName']}: {tar_info[_v_tar]['sourceVariableName']}, QC: {tar_info[_v_tar]['QC']}({qc_thres});{src_info[_v_src]['sourceDataProductName']}: {src_info[_v_src]['sourceVariableName']}"
            gfl_var_info[
                'data_path'] = f"{tar_info[_v_tar]['data_path']}; {src_info[_v_src]['data_path']}"
            gfl_var_info['with_gaps'] = _v_tar

        exp_settings['dataset'][dataset_key]['variables'][
            _v_tar_gfld] = gfl_var_info
        exp_settings['dataset'][dataset_key]['variables'].pop(_v_tar)
    return exp_settings


def write_netcdf(ds, exp_settings, site, resampled=None):
    if resampled is None:
        data_freq = exp_settings["temporal_resolution"]
        ds = add_date_variables(ds)
    else:
        data_freq = resampled
    ncdir = os.path.join(exp_settings['OutPath']['nc_file'], data_freq)
    os.makedirs(ncdir, exist_ok=True)
    ncfile_name = os.path.join(
        ncdir, site + "." + exp_settings['start_date'].split('-')[0] + "." +
        str(int(exp_settings['end_date'].split('-')[0])) + '.' + data_freq +
        '.nc')
    ds.to_netcdf(ncfile_name, mode='w')
    shut.log_and_print_sep(
        f'saved nc data: {data_freq} : {ncfile_name}'.center(150, '-'))
    shut.log_and_print_sep()
    return


def resample_data(_data, _res):

    if _res == 'monthly':
        _data_out = _data.resample(time='1M').mean()
    elif _res == 'annual':
        _data_out = _data.resample(time='1Y').mean()
    else:
        sys.exit(f'{_res} resampling is not implemented')
    return _data_out


def harmonize_data_attrs(site_info, data, _config, resampled=None):
    site = site_info['site_ID']
    start_date = _config["start_date"]
    end_date = _config["end_date"]

    #%Sort variable alphabeticatly
    var_selected_data = sorted(list(data.keys()))
    data = data[var_selected_data]
    if resampled is None:
        data_freq = _config["temporal_resolution"]
    else:
        data_freq = resampled

    #%% Add metadata
    for source_ in _config["sel_datasets"]:
        var_info = _config["dataset"][source_]['variables']
        var_selected_src = list(
            _config["dataset"][source_]['variables'].keys())
        var_selected = list(
            set(var_selected_src).intersection(set(var_selected_data)))

        for var_ in var_selected:
            sourceDataProductName = var_info[var_]['sourceDataProductName']
            if 'time' in data[var_].dims and 'depth_soilGrids' in data[
                    var_].dims:
                nameLong = f"{data_freq} {var_info[var_]['nameLong']} from {sourceDataProductName}"
            else:
                nameLong = f"{var_info[var_]['nameLong']} from {sourceDataProductName}"

            attrs_dict = dict(
                bounds=str([np.nanmin(data[var_]),
                            np.nanmax(data[var_])]),
                units=var_info[var_]['variableUnit'],
                short_name=var_info[var_]['nameShort'],
                long_name=nameLong,
                sourceDataProductName=sourceDataProductName,
                publication=var_info[var_]['publication'],
                partitioning=var_info[var_]['partitioning'],
                sourceVariableName=var_info[var_]['sourceVariableName'],
                description=var_info[var_]['description'],
                data_path=var_info[var_]['data_path'])

            data[var_] = data[var_].assign_attrs(
                collections.OrderedDict(sorted(attrs_dict.items())))
            shut.log_and_print('attrs',
                               site,
                               data[var_].sourceDataProductName,
                               var_,
                               data[var_].units,
                               data_freq,
                               tar_label='variable',
                               src_label='units')

    #%% Transform time dim to days since time dimension
    if resampled is None:
        sdate = _config['start_date'].split('-')
        if _config["temporal_resolution"] == 'daily':
            date_base = datetime(int(sdate[0]), int(sdate[1]), int(sdate[2]),
                                 00, 00, 00)
            date_arr = np.array([
                date_base + timedelta(days=i) for i in range(len(data['time']))
            ])
            data['time'] = date_arr.astype(np.datetime64)
        if _config["temporal_resolution"] == 'hourly':
            date_base = datetime(int(sdate[0]), int(sdate[1]), int(sdate[2]),
                                 00, 00, 00)
            date_arr = np.array([
                date_base + timedelta(hours=i)
                for i in range(len(data['time']))
            ])
            data['time'] = date_arr.astype(np.datetime64)

    #%% Add global attributes
    CurrentScript = os.path.basename(__file__)
    data = data.assign_attrs(
        title=f'site dataset for {_config["exp_name"]}',
        FLUXNET_version=_config["FLUXNET_version"],
        temporal_resolution=data_freq,
        last_disturbance_on=_config['last_disturbance_on'],
        PFT=_config['sitePFT'],
        PFT_extractor=_config['PFTextractor'],
        created_by='Sujan Koirala',
        contact='skoirala@bgc-jena.mpg.de',
        SITE_ID=site_info["site_ID"],
        latitude=site_info["latitude"],
        longitude=site_info["longitude"],
        start_year=_config['start_date'].split('-')[0],
        end_year=_config['end_date'].split('-')[0],
        creation_date=datetime.now().strftime("%d-%m-%Y %H:%M"),
        creation_script=CurrentScript)
    shut.log_and_print_sep()
    return data


def close_logger():
    log = logging.getLogger()  # root logger - Good to get it only once.
    for hdlr in log.handlers[:]:  # remove the existing file handlers
        if isinstance(hdlr, logging.FileHandler):  #fixed two typos here
            log.removeHandler(hdlr)
            hdlr.close()
    return


def compile_site_data(site, exp_settings=None):
    log_path = os.path.join(
        exp_settings['OutPath']['log'],
        f"{site}_{exp_settings['OutPath']['expName']}.log")
    filehandler = logging.FileHandler(log_path, 'w')
    formatter = logging.Formatter(
        '%(levelname)s::%(filename)s::\n\t%(funcName)s::%(lineno)d::\n\t%(message)s'
    )
    filehandler.setFormatter(formatter)
    log = logging.getLogger()
    log.addHandler(filehandler)  # set the new handler
    log.setLevel(logging.INFO)
    # logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
    shut.log_and_print_sep()
    shut.log_and_print_sep(
        f'Initialization of the data processing: {site}'.center(150, ' '))
    shut.log_and_print_sep()
    logging.info(f'Processing data for site {site}')
    #%% Get site information
    site_info = get_site_info(site)

    #%% get all data
    data_all = {}
    sitePFT = ''
    PFTextractor = 'None'
    last_disturbance_on = ''
    for extractor in exp_settings["sel_datasets"]:
        extractor_name = exp_settings["dataset"][extractor]["extractor"]
        logging.info(
            f'Getting dataset| configuration: {extractor}, extractor: {extractor_name}, site: {site_info["site_ID"]}'
        )
        extractor_mod = load_extractor(extractor_name)
        extractor_function = getattr(extractor_mod, 'extract')
        provided_data = extractor_function(site_info=site_info,
                                       config=exp_settings,
                                       dataset=extractor)
        if provided_data is not None:
            data_all[extractor] = provided_data
            if (len(sitePFT) == 0) & ('PFT' in provided_data.attrs):
                sitePFT = f"{provided_data.attrs['PFT']}"
                PFTextractor = f"{extractor} using {extractor_name}"
            if (len(last_disturbance_on) == 0) & ('last_disturbance_on'
                                                  in provided_data.attrs):
                last_disturbance_on = f"{provided_data.attrs['last_disturbance_on']}"
        shut.log_and_print_sep()
    if exp_settings['do_gap_fill'] == True:
        for gfik, gfiv in exp_settings["gap_fill"].items():
            if gfiv['source'] in data_all.keys(
            ) and gfiv['target'] in data_all.keys():
                data_gapfilled = cliff_gapfill(fill_src=gfiv['source'],
                                               fill_tar=gfiv['target'],
                                               site_info=site_info,
                                               data_dict=data_all,
                                               qc_thres=gfiv['qc_thres'],
                                               config=exp_settings).gapfill()
                exp_settings = add_gapfill_settings(exp_settings,
                                                    fill_src=gfiv['source'],
                                                    fill_tar=gfiv['target'],
                                                    dataset_key=f'{gfik}_gfld',
                                                    qc_thres=gfiv['qc_thres'])
                data_all[f'{gfik}_gfld'] = data_gapfilled
            else:
                logging.warning(
                    f"Cannot do gap filling for {gfik}: do_gap_fill is set to {exp_settings['do_gap_fill']} in experiment[json], but one or both of {gfiv['source']} and {gfiv['target']} are not in selected datasets. change sel_datasets, gap filling source or target, or set do_gap_fill to {not exp_settings['do_gap_fill']}."
                )

    # add pft to the settings to write to netcdf

    exp_settings['last_disturbance_on'] = last_disturbance_on
    exp_settings['sitePFT'] = sitePFT
    exp_settings['PFTextractor'] = PFTextractor
    #%% Merge all data
    logging.info(f'Merging all data: {site_info["site_ID"]}')
    try:
        data_values = data_all.values()
        if len(data_values) == 0:
            exit_msg = f"None of the selected datasets provide any data at {exp_settings['temporal_resolution']} resolution. Selected datasets were [{','.join(exp_settings['sel_datasets'])}]"
            logging.warning(exit_msg)
            raise Exception(exit_msg)
        else:
            ds = xr.merge(data_all.values())
    except Exception as e:
        logging.info(
            f"Cannot compile data for: {site_info['site_ID']} for  {site_info['site_ID']}. Error: "
        )
        print(e)
        logging.info(e)
        close_logger()
        return {site: {'status': 'failed'}}

    shut.log_and_print_sep()

    # slice the time series
    ds = ds.sel(
        time=slice(exp_settings["start_date"], exp_settings["end_date"]))

    #%% resample, harmonize metadata, and save
    res_harmon = [None]
    if len(exp_settings["resample_output"]) > 0:
        res_harmon = exp_settings["resample_output"]
        res_harmon.insert(0, None)
    for res in res_harmon:
        if res is not None:
            ds_res = resample_data(ds, res)
        else:
            ds_res = ds
        ds_res_h = harmonize_data_attrs(site_info,
                                        ds_res,
                                        exp_settings,
                                        resampled=res)
        write_netcdf(ds_res_h, exp_settings, site, resampled=res)
        if exp_settings['diagnostic_plots']:
            plot_dianostic_figures(ds_res_h,
                                   site_info,
                                   exp_settings,
                                   site,
                                   resampled=res)
        if res is None:
            res = exp_settings['temporal_resolution']
        shut.log_and_print_sep(
            f'harmonizing, plotting, and saving complete {res}'.center(
                150, '-'))
        shut.log_and_print_sep()

    shut.log_and_print_sep()

    #%% Export site data
    site_dict_d = ds.to_dict(data=False).copy()
    site_dict = {
        site: {
            'start_year': exp_settings['start_date'].split('-')[0],
            'end_year': exp_settings['end_date'].split('-')[0],
            'latitude': site_info['latitude'],
            'longitude': site_info['longitude']
        }
    }
    if 'PFT' in list(site_dict_d['attrs'].keys()):
        site_dict[site]['PFT'] = site_dict['attrs']['PFT']
    # close datasets
    ds_res.close()
    ds.close()

    shut.log_and_print_sep(
        f'Successful data processing for site: {site}'.center(150, ' '))
    shut.log_and_print_sep()

    close_logger()

    return site_dict


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Utilities: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('workflow: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)

# %%

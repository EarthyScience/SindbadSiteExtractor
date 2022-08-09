#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
import pandas as pd
import xarray as xr
import numpy as np

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()

    src_data = []

    src_data_dict = {}

    for tar_name in bxtr.vars_list:


        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        if '_PFT' in tar_name:
            PFT_dat = src_df.loc[src_df['Site_ID'] == bxtr.site].drop(
                ['Site_ID'], axis=1)
            if len(PFT_dat) == 0:
                sel_data = 'undefined'
            else:
                sel_data = PFT_dat['PFT'].values[0]
            print(sel_data)
        else:
            sel_data = np.array(src_df.loc[src_df['Site_ID'] == bxtr.site].drop(
                ['Site_ID'], axis=1),
                                dtype=float)
            data = xr.Dataset({
                tar_name:
                xr.DataArray(data=np.nan, dims=['time'], coords={'time': date_})
            })

            if len(sel_data) > 0:
                sel_data[sel_data == -9999] = np.nan
                sel_data = np.nanmedian(sel_data)
                data[tar_name].values = np.ones_like(
                    data[tar_name].values) * sel_data
            else:
                logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")

            sel_data = bxtr.convert_units(data, tar_name)

        src_data_dict[tar_name] = sel_data


    # merge the data manually (rough method but works.. be careful of hardcoded variable names. things may break if there are changes in extractor's json settings)

    data = xr.where(np.isfinite(src_data_dict['agb_merged']['agb_merged']), src_data_dict['agb_merged']['agb_merged'], src_data_dict['agb_gloBbiomass']['agb_gloBbiomass'])
    tar_name = 'agb_merged'
    # needs renaming and resetting the attributes because xarray likes deleting things
    bxtr.log_var_start(tar_name)
    data = data.rename(tar_name)
    data = bxtr.convert_units(data, tar_name)

    src_data_dict[tar_name][tar_name] = data.copy()


    tar_name = 'agb_merged_PFT'
    # needs renaming and resetting the attributes because xarray likes deleting things
    bxtr.log_var_start(tar_name)
    data = data.rename(tar_name)

    pft = src_data_dict[tar_name]
    if pft in bxtr.vars[tar_name]['PFT_types']:
        data.values = data.values * bxtr.vars[tar_name]['AGB_scalar']
        if np.sum(np.isnan(data.values.flatten())) == len(data.values.flatten()):
            data.values=np.ones_like(data.values) * bxtr.vars[tar_name]['AGB_scalar']
        # data.values = np.ones_like(data.values) * bxtr.vars[tar_name]['AGB_scalar']

    src_data_dict[tar_name] = data

    # only include the merged variables in this extractor to avoid confusion with agb.py. Otherwise, one could just do src_data_dict.values() and include everything
    src_dataset = bxtr.merge_and_format([src_data_dict['agb_merged_PFT'], src_data_dict['agb_merged']])
    return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Provider at: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
import pandas as pd
import xarray as xr
import utils.setup_config as set_conf
import os
import logging
logger = logging.getLogger(__name__)
from extractors.BasexTractor import BasexTractor


def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    in_sub_path = set_conf.get_cliff_dirname(bxtr.version,
                                             bxtr.temporal_resolution)
    date_ = bxtr.get_date_vec()
    src_data = []

    for tar_name in bxtr.vars_list:
        site_data_path = os.path.join(
            bxtr.vars[tar_name]['data_path'], in_sub_path, bxtr.source,
            f'{bxtr.site}.{bxtr.start_date.split("-")[0]}.{bxtr.end_date.split("-")[0]}.{bxtr.temporal_resolution}.csv'
        )
        src_df = pd.read_csv(site_data_path)
        src_name = bxtr.vars[tar_name]['sourceVariableName']
        for ext in ['_DayTime', '_DayMean', '_DayMin', '_DayMax']:
            if ext in tar_name:
                src_name = src_name + ext

        if src_name in src_df:
            sel_data = src_df[src_name].values
            bxtr.log_var_start(tar_name)
            data = xr.Dataset({
                tar_name:
                xr.DataArray(data=sel_data,
                             dims=['time'],
                             coords={'time': date_})
            })
            data = bxtr.convert_units(data, tar_name)
            src_data.append(data)
            bxtr.log_var_end(data, tar_name, None)
        else:
            logger.warning(
                f'{bxtr.site}: target: {tar_name}, src: {src_name} was not found in the csv data frame in {site_data_path}'
            )

    src_dataset = bxtr.merge_and_format(src_data)
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
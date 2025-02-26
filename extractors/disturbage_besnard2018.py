#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A data provicer for forest age and disturbance
Created on Nov 23 2021

@author: dr-ko
"""
import sys, os
sys.path.append(os.path.join(os.getcwd(), '../'))
import pandas as pd
import xarray as xr
import numpy as np

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)

from extractors.BasexTractor import BasexTractor


def extract(dataset, site_info, config):

    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()

    src_data = []

    for tar_name in bxtr.vars_list:
        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        last_disturbance_on = src_df.loc[src_df['Site_ID'] == bxtr.site][
            'Plantation_Date_max'].values.astype(np.datetime64)

        bxtr.log_var_start(tar_name)

        data = xr.Dataset({
            tar_name:
            xr.DataArray(data=0, dims=['time'], coords={'time': date_})
        })

        if len(last_disturbance_on) == 0:
            last_disturbance_on = ['undisturbed']
        else:
            forest_age_data = np.concatenate([
                (day_ - last_disturbance_on).astype('int') / 365.25
                for day_ in date_
            ])
            forest_age_data = xr.DataArray(forest_age_data,
                                           dims=['time'],
                                           coords={'time': date_})
            forest_age_data = forest_age_data.where(forest_age_data == 0)
            forest_age_data = forest_age_data.where(
                np.isfinite(forest_age_data), 1)
            forest_age_data.values = 1 - forest_age_data.values

            data[tar_name].values = forest_age_data.values

        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

    src_dataset = bxtr.merge_and_format(src_data)

    src_dataset = src_dataset.assign_attrs(
        last_disturbance_on=last_disturbance_on[0])

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

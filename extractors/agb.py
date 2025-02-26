#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
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

    for tar_name in bxtr.vars_list:

        bxtr.log_var_start(tar_name)

        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        sel_data = np.array(src_df.loc[src_df['Site_ID'] == bxtr.site].drop(
            ['Site_ID'], axis=1),
                            dtype=float)

        data = xr.Dataset({
            tar_name:
            xr.DataArray(data=np.nan, dims=['time'], coords={'time': date_})
        })

        sel_data[sel_data == -9999] = np.nan
        dval = data[tar_name].values
        if len(sel_data) > 0:
            if '_globBiomass' in tar_name:
                if '_std' in tar_name:
                    sel_data = np.nanstd(sel_data)
                else:
                    sel_data = np.nanmean(sel_data)
                for dI, _date in enumerate(date_):
                    if _date.astype(object).year == 2011:
                        dval[dI-1] = sel_data
                        break
            else:
                for _sel_data in sel_data:
                    yr = _sel_data[0]
                    yr_bio = _sel_data[1]
                    for dI, _date in enumerate(date_):
                        if _date.astype(object).year == yr + 1:
                            dval[dI-1] = yr_bio
                            break
            data[tar_name].values = dval
        else:
            logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")

        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

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
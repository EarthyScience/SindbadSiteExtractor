#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: dr-ko
"""
from fluxcom.core.variables import Variable
from fluxcom.providers import JungQCProvider
from fluxcom.providers import EddyProvider
from fluxcom import transformers

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):

    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    # the following two lines stops the extractor from running when data at a higher resolution than daily is requested. This should be inserted in every extractor to indicate the highest resolution of the data available.

    if not bxtr.is_resolution_supported(provided_reso='daily'):
        return None

    src_prov_j = JungQCProvider(cubepath=bxtr.flx_cubepath,
                              version=bxtr.version,
                              site=bxtr.site)
    src_vars_j = [_var.name for _var in src_prov_j.variables]

    src_prov_ec = EddyProvider(cubepath=bxtr.flx_cubepath,
                                               version=bxtr.version,
                                               site=bxtr.site,
                                               NEE_partitioning_method=None,
                                               transforms=transformers.hourly_to_daily)
    src_vars_ec = [_var.name for _var in src_prov_ec.variables]
    src_data = []
    for tar_name in bxtr.vars_list:
        src_name = bxtr.vars[tar_name]['sourceVariableName']
        if src_name in src_vars_j:
            bxtr.log_var_start(tar_name)

            data_j = src_prov_j.get_data(
                        Variable(src_name,
                         units=bxtr.vars[tar_name]['sourceVariableUnit'],
                         partitioning=bxtr.vars[tar_name]['partitioning']))

            if src_name in src_vars_ec:
                src_name_ec = src_name
            else:
                src_name_ec = 'NEE_QC'
                
            data_ec = src_prov_ec.get_data(
                        Variable(src_name_ec,
                        units=bxtr.vars[tar_name]['sourceVariableUnit'],
                        partitioning=bxtr.vars[tar_name]['partitioning']))

            data = data_j * data_ec
            data = data.rename(tar_name)

            data = bxtr.convert_units(data, tar_name)

            src_data.append(data)

            bxtr.log_var_end(data, tar_name, src_prov_j.transforms)

            src_prov_j.transforms = []

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
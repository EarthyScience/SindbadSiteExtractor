#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
from fluxcom.core.variables import Variable
from fluxcom.providers import JungQCProvider

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):

    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    # the following two lines stops the extractor from running when data at a higher resolution than daily is requested. This should be inserted in every extractor to indicate the highest resolution of the data available.

    if not bxtr.is_resolution_supported(provided_reso='daily'):
        return None

    src_prov = JungQCProvider(cubepath=bxtr.flx_cubepath,
                              version=bxtr.version,
                              site=bxtr.site)
    src_vars = [_var.name for _var in src_prov.variables]
    src_data = []
    for tar_name in list(bxtr.vars.keys()):
        src_name = bxtr.vars[tar_name]['sourceVariableName']
        if src_name in src_vars:
            bxtr.log_var_start(tar_name)

            data = src_prov.get_data(
                Variable(src_name,
                         units=bxtr.vars[tar_name]['sourceVariableUnit'],
                         partitioning=bxtr.vars[tar_name]['partitioning']))
            data = data.rename(tar_name)

            data = bxtr.convert_units(data, tar_name)

            src_data.append(data)

            bxtr.log_var_end(data, tar_name, src_prov.transforms)

            src_prov.transforms = []

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
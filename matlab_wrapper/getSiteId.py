#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 09:29:33 2021

@author: simon
"""
from fluxcom.providers import eddy_covariance as ec
from fluxcom.variable import Variable


def getSiteID(cubepath, version, varname):
    Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = cubepath,
                                                                  version  = version)
    Climate_FLUXNET_data = Climate_FLUXNET_prov.get_data(Variable('P'))
    
    return {'site_id': Climate_FLUXNET_data.site.values}
    
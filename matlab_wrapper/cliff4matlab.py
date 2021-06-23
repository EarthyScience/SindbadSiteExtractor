from fluxcom.providers import eddy_covariance as ec
from fluxcom.variable import Variable
from fluxcom.providers.transformers import hourly_to_daily
import numpy as np

def FLUXCOMdata4cliff(cubepath, version, site, varname, resolution):
    
    ## Import data cube
    if resolution == 'daily':
        Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = cubepath,
                                                                  version  = version, 
                                                                  transforms = hourly_to_daily,                                                     
                                                                  site     = site)
    elif resolution == 'hourly':
            Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = cubepath,
                                                                      version  = version, 
                                                                      site     = site)
    else:
        raise RuntimeError(f'Temporal resolution {resolution} is not implemented')
        
    ## Extract data for each variable
    Climate_FLUXNET_vars = []
    for var_ in varname:
       Climate_FLUXNET_vars.append(Variable(var_))   
    Climate_FLUXNET_data = Climate_FLUXNET_prov.get_data(Climate_FLUXNET_vars)
    
    ## Get latitude and longitude
    lat = Climate_FLUXNET_data.tower_lat.values
    lon = Climate_FLUXNET_data.tower_lon.values

    ## Get julian dates
    date_ = np.arange(np.datetime64('1592-10-15'),  Climate_FLUXNET_data.time[-1:].values[0] + np.timedelta64(1,'D') ,dtype="M8[D]")
    days_since = np.arange(0, len(date_))
    julian_day = days_since[np.in1d(date_, Climate_FLUXNET_data.time)].astype(int)
    year = np.array([str(i).split('-')[0] for i in date_])
    year = year[np.in1d(date_, Climate_FLUXNET_data.time)].astype(int)
    month = np.array([str(i).split('-')[1] for i in date_])
    month = month[np.in1d(date_, Climate_FLUXNET_data.time)].astype(int)
    day = np.array([str(i).split('-')[2] for i in date_])
    day = day[np.in1d(date_, Climate_FLUXNET_data.time)].astype(int)

    ## Get hour
    if resolution == 'daily':
        hour = np.nan
    elif resolution == 'hourly':
        hour = np.tile(Climate_FLUXNET_data.hour +1, Climate_FLUXNET_data.time.shape)
        year  = np.repeat(year, Climate_FLUXNET_data.hour.shape)
        month  = np.repeat(month, Climate_FLUXNET_data.hour.shape)
        day  = np.repeat(day, Climate_FLUXNET_data.hour.shape)
        julian_day  = np.repeat(julian_day, Climate_FLUXNET_data.hour.shape)        
    else:
        raise RuntimeError(f'Temporal resolution {resolution} is not implemented')
    
    ## Get unit for each variable
    vars_dict = []
    unit_vars = {}
    for var_ in list(Climate_FLUXNET_data.keys()):
        vars_dict.append(var_)
        try:
            if resolution == 'daily':
                unit_vars.update({var_: Climate_FLUXNET_data[var_].attrs['units']})
            elif resolution == 'hourly':
                unit_vars.update({var_: Climate_FLUXNET_data[var_].attrs['unit']})
        except :
            unit_vars.update({var_: np.nan})
            
    return {'data':Climate_FLUXNET_data.to_array().values.reshape(len(varname),-1), 'varname': vars_dict, 'units': unit_vars, 'site_id': site, 'latitude': lat, 'longitude': lon, 'julian_day':julian_day, 'year': year, 'month': month, 'day': day, 'hour': hour}
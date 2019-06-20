from numpy import *
import xarray as xr
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.interpolate as spin
import pandas.plotting._converter as pandacnv   # only necessary due to Pandas 0.21.0 bug with Datetime plotting
pandacnv.register()
from datetime import datetime, timedelta
import os
from ecmwfapi import ECMWFDataServer

"""
Saved data files:
    'erai_SH_analysis.nc':             ERA-Interim, 2018-10-10 - 2019-01-01, step=0 (analysis)
                                         times 0:00, 6:00, 12:00, 18:00
                                       grid: 0.75° x 0.75°, area: 53.0°S 180.0°W 90.0°S 180.0°E
    'erai_NH_analysis.nc':             grid: 0.75° x 0.75°, area: 90.0°N 180.0°W 55.0°N 180.0°E
                                 vars: skt - Skin temperature (K) –> (°C)
                                       t2m - Temperature at 2 meters (K) –> (°C)
                                       u10, v10 - U, V wind components at 10 m (m/s)
                                       si10 - 10-m wind speed from 'u10' and 'v10' (m/s)
                                       
    'erai_SH_forecast.nc':             ERA-Interim, 2018-10-10 - 2019-01-01, steps = 6, 12 (forecast)
                                         times 0:00 and 12:00
                                       grid: 0.75° x 0.75°, area: 53.0°S 180.0°W 90.0°S 180.0°E
    'erai_NH_forecast.nc':             grid: 0.75° x 0.75°, area: 90.0°N 180.0°W 55.0°N 180.0°E
                                 vars: tp - Total precipitation (m) –> Precipitation rate (m/s)
                                       sf - Snowfall (m water equivalent) –> Snowfall rate (m/s)
                                       
                                       
Example usage:
    from numpy import *
    import xarray as xr
    import pandas as pd
    import matplotlib.pyplot as plt
    from datetime import datetime, timedelta
    import os
    
    os.chdir('/Users/Ethan/Documents/Research/Git/SnowBlower/')
    sys.path.append('/Users/Ethan/Documents/Research/Git/SnowBlower/')
    import erai_tools as et
    
    erai_SH_forecast = et.load_ecmwf('Data/ERA_Interim_processed/','erai_SH_forecast.nc')
    erai_total_precip_maud_rise = et.create_reanalysis_series(erai_SH_forecast,param_name='tp',
                                                              nearest_to_lat_lon=(-65,0))
                    
                                                              
    # create along-track time series
    data_dir = 'Data/ERA_Interim_processed/'
    erai_analysis_filename = 'erai_SH_analysis.nc'
    erai_forecast_filename = 'erai_SH_forecast.nc'
    
    lats = [-65,-66,-67]
    lons = [5,6,7]
    datetimes = [datetime(2018,12,25,4),datetime(2018,12,25,8),datetime(2018,12,25,12)]
    
    print('Along-track snowfall rate (in m/s) from the past 2 days at each location:')
    print(along_track(lats,lons,datetimes,data_dir,erai_analysis_filename,erai_forecast_filename,
                      temporal='recent_mean',param='snowfall_rate',prior_days=2.0))
    
    print('Current along-track wind speed (in m/s):')
    print(along_track(lats,lons,datetimes,data_dir,erai_analysis_filename,erai_forecast_filename,
                      temporal='current',param='wind_speed_10_m'))

"""

def ecmwf(date_range='1979-01-01/to/2017-08-31',area='-40/-90/-90/90',type='an',step='0',time='00/06/12/18',
          params=['msl','t2m','skt'],output_filename=None):
    """ Submits MARS request to retrieve ERA-Interim reanalysis fields as netCDF file.

    Arguments:
        date_range: for daily fields, format as, e.g., '1979-01-01/to/2017-08-31'
                    for monthly means of daily means, use [datetime(start_yr,start_mo,1),datetime(end_yr,end_mo,1)]
        area: subsetting area, format '-40/-90/-90/90' (N/W/S/E)
        type: 'an' for analysis or 'fc' for forecast
        step: '0' for analysis only, '6/12' or '3/6/9/12' for 6-hourly or 3-hourly forecasts from 0000 and 1200 UTC
              or None for monthly means (regardless, it will be ignored)
        time: analysis times, e.g. '00/06/12/18' for all analyses, or '00/12' if retrieving forecasts only
              or None for monthly means (regardless, it will be ignored)
        params: parameter abbreviations, to be translated into GRIB and Table 2 codes - see below for those available
                note: to find new codes, use parameter database: http://apps.ecmwf.int/codes/grib/param-db/
                      or use web interface and check "View the MARS request"
        output_filename: desired path + filename including '.nc' extension, to save locally
                         or None to save to temporary storage; download from: http://apps.ecmwf.int/webmars/joblist/
                note: if not downloading locally, cancel call using Ctrl-C after "Request is queued" appears
                      (otherwise file will be deleted almost instantly from ECMWF servers)

    None: cancelling call (Ctrl-C) after "Request is queued" appears is fine. It will prevent local download, though.

    Note: private login key required. See documentation for instructions on creating local login key.

    Note: file size limit is probably 20 GB. Check here: https://software.ecmwf.int/wiki/display/WEBAPI/News+feed

    Limited web API access:
        http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/
        http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/

    Documentation:
        https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
        https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
        https://software.ecmwf.int/wiki/display/UDOC/MARS+user+documentation
        https://software.ecmwf.int/wiki/display/UDOC/MARS+keywords
        http://apps.ecmwf.int/codes/grib/param-db

    Reference: Dee et al. 2011

    """
    param_codes = ''
    for param_idx, param in enumerate(params):
        # analysis parameters
        if   param == 't2m':  param_codes += '167.128'  # 2 metre temperature (K)
        elif param == 'sst':  param_codes +=  '34.128'  # Sea surface temperature (K)
        elif param == 'skt':  param_codes += '235.128'  # Skin temperature (K)
        elif param == 'd2m':  param_codes += '168.128'  # 2 metre dewpoint temperature (K)
        elif param == 'msl':  param_codes += '151.128'  # Mean sea level pressure (Pa)
        elif param == 'sp':   param_codes += '134.128'  # Surface pressure (Pa)
        elif param == 'u10':  param_codes += '165.128'  # 10 metre U wind component (m/s)
        elif param == 'v10':  param_codes += '166.128'  # 10 metre V wind component (m/s)
        elif param == 'si10': param_codes += '207.128'  # 10 metre wind speed (m/s) [NOTE: in monthly means only]
        elif param == 'blh':  param_codes += '159.128'  # Boundary layer height (m)
        elif param == 'lcc':  param_codes += '186.128'  # Low cloud cover (fractional coverage, 0 to 1)
        elif param == 'tcc':  param_codes += '164.128'  # Total cloud cover (fractional coverage, 0 to 1)
        elif param == 'rsn':  param_codes +=  '33.128'  # Snow density in snow layer (kg/m^3)
        elif param == 'sd':   param_codes += '141.128'  # Snow depth in snow layer (m of water equivalent)
        elif param == 'sr':   param_codes += '173.128'  # Climatological aerodynamic land surface roughness length (m)
        elif param == 'tsn':  param_codes += '238.128'  # Temperature of snow layer (K)
        # forecast parameters (* indicates accumulated field; note downward fluxes are positive)
        elif param == 'sf':   param_codes += '144.128'  # Snowfall (m of water equivalent) *
        elif param == 'sshf': param_codes += '146.128'  # Surface sensible heat flux (J/m^2) *
        elif param == 'slhf': param_codes += '147.128'  # Surface latent heat flux (J/m^2) *
        elif param == 'ssr':  param_codes += '176.128'  # Surface net solar radiation [shortwave] (J/m^2) *
        elif param == 'str':  param_codes += '177.128'  # Surface net thermal radiation [longwave] (J/m^2) *
        elif param == 'strd': param_codes += '175.128'  # Surface thermal radiation [longwave] downwards (J/m^2) *
        elif param == 'e':    param_codes += '182.128'  # Evaporation (m of water equivalent) *
        elif param == 'tp':   param_codes += '228.128'  # Total precipitation (m) *
        elif param == 'iews': param_codes += '229.128'  # Instantaneous eastward turbulent surface stress (N/m^2)
        elif param == 'inss': param_codes += '230.128'  # Instantaneous northward turbulent surface stress (N/m^2)
        if param_idx < len(params)-1: param_codes += '/'

    retrieve_dict = {
        "class":"ei",
        "dataset":"interim",
        "expver":"1",
        "format":"netcdf",
        "grid":"0.75/0.75",
        "levtype":"sfc",
        "param":param_codes,
        "type":type,
        'area':area,
        "target":output_filename,
        "use":'frequent',
    }

    # monthly means of daily means
    if len(date_range) == 2:
        retrieve_dict['stream'] = 'moda'
        final_date_range = ''
        working_month = date_range[0]
        while working_month < date_range[1]:
            final_date_range += working_month.strftime('%Y%m%d')
            final_date_range += '/'
            working_month += relativedelta(months=+1)
        final_date_range += date_range[1].strftime('%Y%m%d')
        retrieve_dict['date'] = final_date_range

    # daily fields
    else:
        retrieve_dict['stream'] = 'oper'
        retrieve_dict['date'] = date_range
        retrieve_dict['step'] = step
        retrieve_dict['time'] = time

    server = ECMWFDataServer()
    server.retrieve(retrieve_dict)


def load_ecmwf(data_dir,filename,datetime_range=None,lat_range=None,lon_range=None,
               export_to_dir=None,export_filename=None,export_chunks=True,verbose=False,super_verbose=False):
    """ Opens ERA-Interim reanalysis data files downloaded in netCDF format with a custom grid.

    Secondary use:
        Run this routine on newly downloaded files to calculate derived quantities and de-accumulate forecasts. Use
        argument <<export_to_dir>> to export new version, then manually delete original.

    Args:
        data_dir: directory of data file
        filename: filename including extension
        datetime_range: None or [Datetime0,Datetime1] or [Datestring0,Datestring1] to subset fields
            note: to slice with open right end, e.g., use [Datetime0,None]
            note: selection is generous, so ['2016-1-1','2016-1-1'] will include all hours on January 1, 2016
            note: example of Datestring: '2016-1-1-h12'
        lat_range: None or [lat_N,lat_S] to subset fields (!!! - descending order - !!!)
        lon_range: None or [lon_W,lon_E] to subset fields
        export_to_dir: None or directory to export new netCDF containing derived quantities, modified variables
            note: in this case, will not return Dataset, and will ignore datetime_range when export_chunks is True
            note: this must be a *different* directory than data_dir!
        export_filename: None or new filename to use when exporting (including extension)
        export_chunks: True or False (True for first call on a large file; False when calling on chunks of that file)
        verbose: True or False
        super_verbose: True or False (print every processing time step)

    Returns:
        all_data: xarray Dataset with coordinates (time,lats,lons); examples of accessing/slicing follow:
            all_data.loc[dict(time='2016-1-1')]                            to extract without slicing
            all_data.sel(lats=slice(-60,-70))                              to slice all variables
            all_data['skt'].values                                         to convert to eager NumPy array
            all_data['skt'][0,:,:]                                         to slice data using indices (t=0)
            all_data['skt'].loc['2016-1-1':'2016-2-1',-60:-70,0:10]        to slice data using values
            all_data['skt']['latitude']                                    to get view of 1-D coordinate
            all_data['skt']['time']                                        NumPy Datetime coordinate
            all_data['skt']['doy']                                         fractional day-of-year coordinate
            pd.to_datetime(all_data['skt']['time'].values)                 useable Datetime version of the above
            all_data['skt'].attrs['units']
            all_data['skt'].attrs['long_name']

    Note: as shown above, 'doy' (fractional day-of-year) is included as a secondary coordinate with dimension 'time'.

    The following derived quantities are calculated here:
        'q2m': 2-m specific humidity from 'msl' and 'd2m'
        'si10': 10-m wind speed from 'u10' and 'v10' (evaluated lazily using Dask only if export_to_dir is None)

    """

    # export mode may require splitting numerical processing into chunks
    max_chunk = 0.5  # in GB, maximum file size to attempt to process in memory
    if export_to_dir is not None and export_chunks:
        file_size = os.path.getsize(data_dir + filename)/10e8  # file size in GB
        if file_size > max_chunk:
            num_chunks = int(ceil(file_size/max_chunk))
            all_data = xr.open_dataset(data_dir + filename)
            num_times = all_data.dims['time']
            times_per_chunk = int(ceil(num_times/num_chunks))
            all_times = all_data['time'].values
            all_data.close()

            # process and save data in chunks
            slice_start_indices = arange(0,num_times,times_per_chunk)
            just_filename,just_extension = os.path.splitext(filename)
            for chunk_counter, start_idx in enumerate(slice_start_indices):
                end_idx = start_idx + times_per_chunk - 1
                if end_idx >= len(all_times): end_idx = -1  # for final chunk, use last time index
                dt_range = [str(all_times[start_idx]),str(all_times[end_idx])]
                if verbose: print('>> Processing chunk {0} of {1} from {2} to {3}'
                                  ''.format(chunk_counter+1,len(slice_start_indices),*dt_range))
                load_ecmwf(data_dir,filename,datetime_range=dt_range,lat_range=lat_range,lon_range=lon_range,
                           export_filename='{0}_chunk_{1:03d}{2}'.format(just_filename,chunk_counter+1,just_extension),
                           export_to_dir=export_to_dir,export_chunks=False,verbose=verbose)

            # open all chunks and concatenate as Dataset
            if verbose: print('>> Opening all chunks of {0}'.format(filename))
            all_data = xr.open_mfdataset(export_to_dir + '{0}_chunk_*{1}'.format(just_filename,just_extension),
                                         concat_dim='time',chunks={'time':100})
            bypass_normal_open = True
        else:
            bypass_normal_open = False
    else:
        bypass_normal_open = False

    if not bypass_normal_open:
        if verbose: print('>> Opening {0}'.format(filename))
        all_data = xr.open_dataset(data_dir+filename,chunks={'time':100})   # O(100 MB) per chunk

    if 'longitude' in all_data and 'latitude' in all_data:
        all_data = all_data.rename({'latitude':'lats','longitude':'lons'})

    if datetime_range is not None:
        all_data = all_data.sel(time=slice(datetime_range[0],datetime_range[1]))
    if lat_range is not None:
        all_data = all_data.sel(lats=slice(lat_range[0],lat_range[1]))
    if lon_range is not None:
        all_data = all_data.sel(lons=slice(lon_range[0],lon_range[1]))

    for var in all_data.data_vars:
        if verbose: print('>>>> Examining variable {0}'.format(var))

        if all_data[var].attrs['units'] == 'Pa':
            orig_name = all_data[var].attrs['long_name']
            all_data[var] /= 100.0
            all_data[var].attrs = {'units':'hPa','long_name':orig_name}
        elif all_data[var].attrs['units'] == 'K' and var != 'd2m':
            orig_name = all_data[var].attrs['long_name']
            all_data[var] -= 273.15
            all_data[var].attrs = {'units':'°C','long_name':orig_name}

        # de-accumulate forecast fields (hours 0 and 12), if not already
        if var in ['tp','e','sf','sshf','slhf','ssr','str','strd'] and 'deaccumulated' not in all_data[var].attrs:
            orig_name = all_data[var].attrs['long_name']
            orig_units = all_data[var].attrs['units']
            time_index = pd.to_datetime(all_data[var].time.values)

            if time_index[0].hour == 0 or time_index[0].hour == 12:
                all_data[var][dict(time=0)] /= 2
                first_step = 1
            else:
                first_step = 0
            if time_index[-1].hour in [3,6,9,15,18,21]:  # handles 3-hourly and 6-hourly steps
                last_step = len(time_index) - 1
            else:
                last_step = len(time_index)

            all_data[var].load() # load Dask array into memory (which means reasonably small chunks are necessary!)
            all_data[var][first_step+1:last_step:2] -= all_data[var][first_step:last_step:2].values

            if   (all_data['time'].values[1] - all_data['time'].values[0]) == timedelta64(6,'h'):
                step_hours = 6
            elif (all_data['time'].values[1] - all_data['time'].values[0]) == timedelta64(3,'h'):
                step_hours = 3
            else:
                raise ValueError('Error from ldp.load_ecmwf(): Forecast time interval is not 3 or 6 hours.')

            seconds_in_step = step_hours * 60 * 60
            all_data[var] /= seconds_in_step

            if var == 'e': all_data[var] *= -1

            all_data[var].attrs['long_name'] = orig_name
            all_data[var].attrs['units'] = orig_units

            if   all_data[var].attrs['units'] == 'm':       all_data[var].attrs['units'] = 'm/s'
            elif all_data[var].attrs['units'] == 'J m**-2': all_data[var].attrs['units'] = 'W/m^2'

            all_data[var].attrs['deaccumulated'] = 'True'

    # calculate 2-m specific humidity from surface pressure and dewpoint temperature, if available
    # uses Equations 7.4 and 7.5 on p. 92 of ECMWF IFS Documentation, Ch. 7:
    #   https://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf
    if 'q2m' not in all_data and 'd2m' in all_data and 'msl' in all_data:
        if verbose: print('>>>> Calculating 2-m specific humidity')

        # constants for Teten's formula for saturation water vapor pressure over water [not ice] (Eq. 7.5)
        # origin: Buck (1981)
        a1 = 611.21 # Pa
        a3 = 17.502 # unitless
        a4 = 32.19  # K
        T_0 = 273.16 # K

        # saturation water vapor pressure; units: Pa
        e_sat_at_Td = a1 * exp(a3 * (all_data['d2m'] - T_0) / (all_data['d2m'] - a4))

        # saturation specific humidity at dewpoint temperature (Eq. 7.4)
        # note conversion of surface pressure from hPa back to Pa
        R_dry_over_R_vap = 0.621981  # gas constant for dry air over gas constant for water vapor, p. 110
        q_sat_at_Td = R_dry_over_R_vap * e_sat_at_Td / (100*all_data['msl'] - (e_sat_at_Td*(1.0 - R_dry_over_R_vap)))

        all_data['q2m'] = q_sat_at_Td
        all_data['q2m'].attrs['units'] = 'kg/kg'
        all_data['q2m'].attrs['long_name'] = 'Specific humidity at 2 m'

    # calculate 10-m wind speed from u, v
    # note: this evaluates lazily using Dask, so expect processing hangs upon computation (instead of load)
    # note: included only if not exporting to a new netCDF file (don't want to take up unnecessary space)
    if 'si10' not in all_data and 'u10' in all_data and 'v10' in all_data and export_to_dir is None:
        if verbose: print('>>>> Calculating 10-m wind speed')
        all_data['si10'] = (all_data['u10']**2 + all_data['v10']**2)**0.5
        all_data['si10'].attrs['units'] = 'm/s'
        all_data['si10'].attrs['long_name'] = '10 metre wind speed'

    # add day-of-year as a secondary coordinate with dimension 'time'
    if 'doy' not in all_data.coords:
        datetime_index = pd.to_datetime(all_data['time'].values)
        doy_index = datetime_index.dayofyear + datetime_index.hour / 24. + datetime_index.minute / 60.
        all_data.coords['doy'] = ('time',doy_index)

    if export_to_dir is not None:
        # set encoding only if exporting to a new netCDF file here (!)
        # remember to do this if exporting to a new netCDF file elsewhere...
        #
        # changing encoding (scale factor and offset) is necessary because the original netCDF file's encoding
        #   results in truncation/loss of precision when applied to the processed variables here (some of which
        #   where divided by large numbers, for instance)
        # these formulae for optimal scale factors and offsets are from:
        #   http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors/
        for var in all_data.data_vars:
            n_bits = 16  # because int16
            var_max = asscalar(all_data[var].max().values)  # .values necessary because of lazy Dask evaluation
            var_min = asscalar(all_data[var].min().values)
            all_data[var].encoding['dtype'] = 'int16'
            all_data[var].encoding['scale_factor'] = (var_max - var_min) / ((2**n_bits) - 1)
            all_data[var].encoding['add_offset'] = var_min + (2**(n_bits - 1) * all_data[var].encoding['scale_factor'])
            all_data[var].encoding['_FillValue'] = -9999

        if export_filename is None: new_filename = filename
        else:                       new_filename = export_filename
        all_data.to_netcdf(export_to_dir + new_filename)
        all_data.close()
    else:
        return all_data


def load_ecmwf_mask(data_dir,filename,var_name='lsm'):
    """ Returns xarray DataArray of mask (e.g. land-sea mask) for ECMWF reanalysis grid (e.g. ERA-Interim).

    Downloaded manually:
        0.75x0.75° ERA-Interim land-sea mask, netCDF: http://apps.ecmwf.int/datasets/data/interim-full-invariant/
            variable name is 'lsm' (0 = sea, 1 = land)

    """
    mask = xr.open_dataset(data_dir + filename)
    mask = mask[var_name].isel(time=0)

    mask = mask.rename({'latitude':'lats','longitude':'lons'})
    lons = mask['lons'].values
    lons[lons > 180.0] -= 360
    mask['lons'] = lons

    return mask


def create_reanalysis_series(dataset,param_name='msl',avg_box=None,nearest_to_lat_lon=None,
                             mask_land=None,mask_sea=None,avg_box_north_of_antarctic_coast=False,circumant_lats=None,
                             min_not_mean=False,max_not_mean=False,abs_val=False,
                             create_climo=False,create_climo_iqr=False,climo_years=None,make_year=None,
                             climo_abs_val=['u10','v10','iews','inss']):
    """ Create a record/index of any reanalysis parameter (<<param_name>>) averaged within a given box (<<avg_box>>).

    Args:
        dataset: xarray Dataset produced by load_ecmwf()
        param_name: abbreviation (key) for parameter of interest
        avg_box: list of length 4, [lon_E,lon_W,lat_S,lat_N], representing extent of averaging box
        nearest_to_lat_lon: None or (lat,lon) tuple to only analyze nearest grid cell, instead of using avg_box
        mask_land: None or xarray DataArray with corresponding land-sea mask (0 = sea, 1 = land) to ignore land
        mask_sea: None or [see above] to ignore sea
        avg_box_north_of_antarctic_coast: if True, mask south of the Antarctic coast
            note: this is an earlier, more crude version of mask_land above
        circumant_lats: vector of Antarctic coast latitudes, corresponding to longitudes of -180 to 179 (inclusive)
                        with a spacing of 1.0 (must be supplied if <<avg_box_north_of_antarctic_coast>> is True)
        min_not_mean: return minimum value found within <<avg_box>>, not the mean
        max_not_mean: return maximum value found within <<avg_box>>, not the mean
        abs_val: use absolute value of data when finding max, min, or mean in avg_box
        create_climo: if True, return only climatology (mean, std) instead of return args below
                    note: climo series format is Pandas Series with fractional DOY index (e.g. 1.0 to 365.75)
        create_climo_iqr: if True, return climatology (median, iqr_25, iqr_75) instead of return args below
                    note: climo series format is Pandas Series with fractional DOY index (e.g. 1.0 to 365.75)
        climo_years: None (to use all available years)
                     or years as [start,end] (inclusive) to specify which to include in climatology
        make_year: None or specific year (integer YYYY) to export Series (either normal indices, or climo) with
                           Datetime index starting from Jan. 1 of that year
                           note: requires one-year-long series; otherwise non-unique datetimes will be generated
        climo_abs_val: list of parameters for which to calculate climatology on absolute value of series

    Returns:
        index: Pandas Series of calculated values with Pandas Timestamp index
               note: Pandas Timestamps are essentially interchangeable with Datetime, but if conversion needed, see:
                     https://pandas.pydata.org/pandas-docs/stable/generated/pandas.Timestamp.html

    """
    if nearest_to_lat_lon is not None:
        index = dataset[param_name].sel(lats=nearest_to_lat_lon[0],lons=nearest_to_lat_lon[1],
                                        method='nearest')
    else:
        if mask_land is not None or mask_sea is not None:
            if mask_land is not None: geo_mask = mask_land; ignore_val = 1.0
            else:                     geo_mask = mask_sea;  ignore_val = 0.0
            data = dataset[param_name].where(geo_mask != ignore_val)
        elif avg_box_north_of_antarctic_coast:
            lon_grid, lat_grid = meshgrid(dataset[param_name]['lons'],dataset[param_name]['lats'])
            geo_mask = lat_grid > circumant_lats[(floor(lon_grid) + 180).astype(int)]
            data = dataset[param_name].where(geo_mask)
        else:
            data = dataset[param_name]

        if abs_val: data = xr.ufuncs.fabs(data)  # maintains lazy evaluation of Dask array

        data = data.sel(lons=slice(avg_box[0],avg_box[1]),lats=slice(avg_box[3],avg_box[2]))

        # note: Dask array computation/conversion triggered here by .compute()
        #       could also keep as Dask array using .persist(), but not sure if this has any advantages
        if min_not_mean:   index = data.min(dim=['lats','lons'],keep_attrs=True,skipna=True).compute()
        elif max_not_mean: index = data.max(dim=['lats','lons'],keep_attrs=True,skipna=True).compute()
        else:              index = data.mean(dim=['lats','lons'],keep_attrs=True,skipna=True).compute()

    if create_climo or create_climo_iqr:
        if param_name in climo_abs_val: index = abs(index)
        if climo_years is None: index_trimmed = index
        else:                   index_trimmed = index.loc[str(climo_years[0]):str(climo_years[1])]
        if create_climo:
            climo = index_trimmed.groupby('doy').mean(dim='time')
            climo_std = index_trimmed.groupby('doy').std(dim='time')
            climo_series = climo.to_pandas()
            climo_std_series = climo_std.to_pandas()
        elif create_climo_iqr:
            climo = index_trimmed.groupby('doy').median(dim='time')
            climo_iqr_25 = index_trimmed.groupby('doy').reduce(stats.iqr,rng=(25,50))
            climo_iqr_75 = index_trimmed.groupby('doy').reduce(stats.iqr,rng=(50,75))
            climo_series = climo.to_pandas()
            climo_iqr_25_series = climo_iqr_25.to_pandas()
            climo_iqr_75_series = climo_iqr_75.to_pandas()
        if make_year is not None:
            ref_datetime = datetime(make_year-1,12,31)
            new_index = array([timedelta(days=int(doy)) + timedelta(hours=int((doy-floor(doy))*24))
                               for doy in climo_series.index]) + ref_datetime
            climo_series.index = new_index
            if create_climo:
                climo_std_series.index = new_index
            elif create_climo_iqr:
                climo_iqr_25_series.index = new_index
                climo_iqr_75_series.index = new_index
        if create_climo:
            return climo_series, climo_std_series
        elif create_climo_iqr:
            return climo_series, climo_iqr_25_series, climo_iqr_75_series

    # convert from xarray to Pandas (because rolling() operations buggy in xarray)
    index = index.to_pandas()
    if make_year is not None:
        new_index = array([dt.replace(year=make_year) for dt in index.index])
        index.index = new_index

    return index


def datetime_to_datenum(datetime_vector):
    return mdates.datestr2num([str(dt) for dt in datetime_vector])


def series_interp(orig_series,new_index,day_interval=(1 / 24)):
    # note: new_index can be a single Datetime or a range of Datetimes, e.g. [start,end]
    interpolator = spin.interp1d(datetime_to_datenum(orig_series.index),orig_series.values,
                                 bounds_error=False,fill_value=NaN)
    if not isinstance(new_index,datetime):
        new_index = arange(*new_index,timedelta(days=day_interval))
    else:
        new_index = [new_index]
    interp_data = interpolator(datetime_to_datenum(new_index))
    return pd.Series(index=new_index,data=interp_data)


def loc_history(lat,lon,this_dt,erai_analysis,erai_forecast,prior_days=2.0):
    erai_analysis_subset = erai_analysis.sel(time=slice(this_dt - timedelta(days=prior_days + 1.0),
                                                         this_dt + timedelta(days=1.0)))
    erai_forecast_subset = erai_forecast.sel(time=slice(this_dt - timedelta(days=prior_days + 1.0),
                                                         this_dt + timedelta(days=1.0)))

    time_range = [this_dt - timedelta(days=prior_days),this_dt]

    recent_sf = create_reanalysis_series(erai_forecast_subset,param_name='sf',nearest_to_lat_lon=(lat,lon))
    recent_ws = create_reanalysis_series(erai_analysis_subset,param_name='si10',nearest_to_lat_lon=(lat,lon))
    recent_t2m = create_reanalysis_series(erai_analysis_subset,param_name='t2m',nearest_to_lat_lon=(lat,lon))
    recent_sf_series = series_interp(recent_sf,time_range)
    recent_ws_series = series_interp(recent_ws,time_range)
    recent_t2m_series = series_interp(recent_t2m,time_range)

    data_dict = {'current':{'snowfall_rate':recent_sf_series[-1],
                            'wind_speed_10_m':recent_ws_series[-1],
                            'temp_2_m':recent_t2m_series[-1]},
                 'recent_mean':{'snowfall_rate':mean(recent_sf_series),
                                'wind_speed_10_m':mean(recent_ws_series),
                                'temp_2_m':mean(recent_t2m_series)},
                 'recent_max':{'snowfall_rate':max(recent_sf_series),
                               'wind_speed_10_m':max(recent_ws_series),
                               'temp_2_m':max(recent_t2m_series)}}
    return data_dict


def along_track(lats,lons,datetimes,data_dir,erai_analysis_filename,erai_forecast_filename,prior_days=2.0,
                temporal='recent_mean',param='snowfall_rate'):
    """
    Args:
        lats, lons: NumPy arrays of along-track latitudes and longitudes
        datetimes: NumPy arrays of along-track Datetime objects
        prior_days: number of days over which to calculate averages or maxima
        temporal: {'current','recent_mean','recent_max'}
        param: {'snowfall_rate','wind_speed_10_m','temp_2_m'}

    """
    erai_analysis = load_ecmwf(data_dir,erai_analysis_filename)
    erai_forecast = load_ecmwf(data_dir,erai_forecast_filename)

    return array([loc_history(lats[t_idx],lons[t_idx],datetimes[t_idx],erai_analysis,erai_forecast,
                              prior_days=prior_days)[temporal][param] for t_idx in range(0,len(datetimes))])


if __name__== "__main__":
    # control flow
    download_ecmwf = True
    process_ecmwf = False
    load_erai_land_sea_mask = True

    # filepaths
    main_dir = os.getcwd()
    era_old_dir = main_dir + '/Data/ERA_Interim_unprocessed/'
    era_new_dir = main_dir + '/Data/ERA_Interim_processed/'

    # submit MARS request for ERA-Interim reanalysis fields
    # note: submit one at a time; cancel using Ctrl-C (or "stop" button) immediately after seeing "Request is queued"
    #       then download using Chrome from: http://apps.ecmwf.int/webmars/joblist/
    #       and save using filenames in comments in folder 'ERA_Interim_unprocessed'
    #       then run 'process_ecmwf' routine below
    if download_ecmwf:
        which_to_download = 3  # change to submit one at a time (see note above) - recommend order 3, 4, 1, 2

        if which_to_download == 1:      # analysis; save as 'erai_SH_analysis.nc'
            ecmwf(date_range='2018-10-10/to/2019-01-01',area='-53/-180/-90/180',output_filename=None,type='an',
                  step='0',time='00/06/12/18',params=['skt','t2m','u10','v10','blh','lcc','tcc','rsn','sd','sr','tsn'])
        elif which_to_download == 2:    # forecast; save as 'erai_SH_forecast.nc'
            ecmwf(date_range='2018-10-10/to/2019-01-01',area='-53/-180/-90/180',output_filename=None,type='fc',
                  step='6/12',time='00/12',params=['tp','sf'])
        elif which_to_download == 3:    # analysis; save as 'erai_NH_analysis.nc'
            ecmwf(date_range='2018-10-10/to/2019-01-01',area='90/-180/55/180',output_filename=None,type='an',
                  step='0',time='00/06/12/18',params=['skt','t2m','u10','v10','blh','lcc','tcc','rsn','sd','sr','tsn'])
        elif which_to_download == 4:    # forecast; save as 'erai_NH_forecast.nc'
            ecmwf(date_range='2018-10-10/to/2019-01-01',area='90/-180/55/180',output_filename=None,type='fc',
                  step='6/12',time='00/12',params=['tp','sf'])

    # process newly downloaded ECMWF reanalysis files (calculate derived quantities, de-accumulate, and re-export)
    # note: once finished, manually delete unprocessed files and any processed chunks
    if process_ecmwf:
        for filename in os.listdir(path=era_old_dir):
            if filename == '.DS_Store': continue
            load_ecmwf(era_old_dir,filename,export_to_dir=era_new_dir,verbose=True)

    # load ERA-Interim land-sea mask
    if load_erai_land_sea_mask:
        erai_mask = load_ecmwf_mask(era_new_dir,'erai_land_sea_mask.nc')



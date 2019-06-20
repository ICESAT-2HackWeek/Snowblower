import numpy as np
import h5py
from astropy.time import Time
import os

def gps2dyr(time):
    """ Convert GPS time to decimal years. """
    return Time(time, format='gps').decimalyear

def track_type(time,lat,tmax=1):
    """
    Separate tracks into ascending and descending.

    Defines tracks as segments with time breaks > tmax,
    and tests whether lat increases or decreases w/time.
    """
    tracks = np.zeros(lat.shape)  # generate track segment
    tracks[0:np.argmax(np.abs(lat))] = 1  # set values for segment
    i_asc = np.zeros(tracks.shape,dtype=bool)  # output index array

    # Loop trough individual segments
    for track in np.unique(tracks):

        i_track, = np.where(track == tracks)  # get all pts from seg

        if len(i_track) < 2: continue

        # Test if lat increases (asc) or decreases (des) w/ time
        i_min = time[i_track].argmin()
        i_max = time[i_track].argmax()
        lat_diff = lat[i_track][i_max] - lat[i_track][i_min]

        # Determine track type
        if lat_diff > 0:  i_asc[i_track] = True

    return i_asc,np.invert(i_asc)  # index vectors

def process(fname,bbox=None):
    """
    Read 1 ATL06 file and output 6 reduced files.

    Extract variables of interest and separate the ATL06 file
    into each beam (ground track) and ascending/descending orbits.
    """

    # Each beam is a group
    group = ['/gt1l','/gt1r','/gt2l','/gt2r','/gt3l','/gt3r']

    # Loop trough beams
    for k,g in enumerate(group):  # 6 groups for 6 beams

        # -----------------------------------#
        # 1) Read in data for a single beam #
        # -----------------------------------#

        # Load variables into memory (more can be added!)
        with h5py.File(fname,'r') as fi:
            lat = fi[g + '/land_ice_segments/latitude'][:]
            lon = fi[g + '/land_ice_segments/longitude'][:]
            h_li = fi[g + '/land_ice_segments/h_li'][:]
            h_li_sigma = fi[g + '/land_ice_segments/h_li_sigma'][:]
            t_dt = fi[g + '/land_ice_segments/delta_time'][:]
            q_flag = fi[g + '/land_ice_segments/atl06_quality_summary'][:]
            dem_h = fi[g + '/land_ice_segments/dem/dem_h'][:]
            bsnow_h = fi[g + '/land_ice_segments/geophysical/bsnow_h'][:]
            bsnow_conf = fi[g + '/land_ice_segments/geophysical/bsnow_conf'][:]
            bsnow_od = fi[g + '/land_ice_segments/geophysical/bsnow_od'][:]
            cloud_flg_asr = fi[g + '/land_ice_segments/geophysical/cloud_flg_asr'][:]
            cloud_flg_atm = fi[g + '/land_ice_segments/geophysical/cloud_flg_atm'][:]
            msw_flag = fi[g + '/land_ice_segments/geophysical/msw_flag'][:]
            t_ref = fi['/ancillary_data/atlas_sdp_gps_epoch'][:]
            orb = np.full_like(h_li,k)

        # ---------------------------------------------#
        # 2) Filter data according region and quality #
        # ---------------------------------------------#

        # Select a region of interest
        if bbox:
            lonmin,lonmax,latmin,latmax = bbox
            bbox_mask = (lon >= lonmin) & (lon <= lonmax) & \
                        (lat >= latmin) & (lat <= latmax)
        else:
            bbox_mask = np.ones_like(lat,dtype=bool)  # get all

        # Only keep good data, and data inside bbox
        mask = (q_flag == 0) & (np.abs(h_li) < 10e3) & (bbox_mask == 1)

        # Update variables
        lat,lon,h_li,h_li_sigma,t_dt,q_flag,dem_h,bsnow_conf,bsnow_h,bsnow_od,cloud_flg_asr,cloud_flg_atm,msw_flag = \
            lat[mask],lon[mask],h_li[mask],h_li_sigma[mask],t_dt[mask],q_flag[mask],dem_h[mask], \
            bsnow_conf[mask],bsnow_h[mask],bsnow_od[mask], \
            cloud_flg_asr[mask],cloud_flg_atm[mask],msw_flag[mask]

        # Test for no data
        if len(h_li) == 0:
            print('ERROR from read_atl06(): granule found with zero length <<h_li>>\n')
            continue

        # -------------------------------------#
        # 3) Convert time and separate tracks #
        # -------------------------------------#

        # Time in GPS seconds (secs sinde 1980...)
        t_gps = t_ref + t_dt

        # Time in decimal years
        t_year = gps2dyr(t_gps)

        # Determine orbit type
        i_asc,i_des = track_type(t_year,lat)

        # -----------------------#
        # 4) Save selected data #
        # -----------------------#

        # Define output file name
        ofile = fname.replace('.h5','_') + 'reduced.h5'

        # Save variables
        with h5py.File(ofile,'a') as f:
            f[g[1:] + '/orbit'] = orb
            f[g[1:] + '/lon'] = lon
            f[g[1:] + '/lat'] = lat
            f[g[1:] + '/t_year'] = t_year
            f[g[1:] + '/t_sec'] = t_gps
            f[g[1:] + '/q_flag'] = q_flag
            f[g[1:] + '/h_li'] = h_li
            f[g[1:] + '/h_li_sigma'] = h_li_sigma
            f[g[1:] + '/dem_h'] = dem_h
            f[g[1:] + '/bsnow_conf'] = bsnow_conf
            f[g[1:] + '/bsnow_h'] = bsnow_h
            f[g[1:] + '/bsnow_od'] = bsnow_od
            f[g[1:] + '/cloud_flg_asr'] = cloud_flg_asr
            f[g[1:] + '/cloud_flg_atm'] = cloud_flg_atm
            f[g[1:] + '/msw_flag'] = msw_flag
            f[g[1:] + '/trk_type'] = i_asc

            print('out ->',ofile)

def process_all(dir):
    """ Processes all .h5 files in the given directory using process().
    """
    all_files = os.listdir(dir)
    if '.DS_Store' in all_files: all_files.remove('.DS_Store')
    remove_these = []
    for fn in all_files:
        if '.h5' not in fn:
            remove_these.append(fn)
    for rfn in remove_these:
        all_files.remove(rfn)
    for fn in all_files:
        process(fn)

def read_h5(fname,vnames=[]):
    """ Simple HDF5 reader. """
    with h5py.File(fname,'r') as f:
        beam_keys = f.keys()
        data = dict()
        for bk in beam_keys:
            data[bk] = dict()
            for vn in vnames:
                data[bk][vn] = f[bk][vn][:]
    return data

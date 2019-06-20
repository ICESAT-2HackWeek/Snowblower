import numpy as np
import h5py
import pandas as pd
import os

def downsample(h5path):
	lon_range = np.arange(-180, 179.9, 0.75)
	lat_range = np.arange(-69.75, -90.1, -0.75)+0.375
	laser_strs = ['gt1r','gt2r','gt3r', 'gt1l','gt2l','gt3l']
	f = h5py.File(h5path,'r')

	df = pd.DataFrame(columns = ['laser', 'lon', 'lat', 'bs_conf', 'bs_conf_x2', 'bs_od', 'bs_od_x2', 'bs_h', 'bs_h_x2',
								 'N','n_bsnow','%_seg', 'h_li', 'mean_timestamp_s', 'mean_timestamp_yr'])

	l=0
	for k in range(len(laser_strs)):
		lastr = laser_strs[k]
		lat = f[lastr+'/lat'][:]
		lon = f[lastr+'/lon'][:]
		bsnow_conf = f[lastr+'/bsnow_conf'][:]
		bsnow_od  = f[lastr+'/bsnow_od'][:]
		bsnow_h = f[lastr+'/bsnow_h'][:]
		dt = f[lastr+'/t_sec'][:]
		dy = f[lastr+'/t_year'][:]
		h_li = f[lastr+'/h_li'][:]        # NOTE: getting questionable NaNs - check this!

		lon[lon>179.625]-=360
		lon_inds = np.digitize(lon, lon_range)
		lat_inds = np.digitize(lat, lat_range)

		for i in range(len(lon_range)):
			for j in range(len(lat_range)):

				# i = 257
				# j = 15

				grab_inds = np.argwhere((lat_inds == j)*(lon_inds==i))
				if len(grab_inds)>0:
					N = len(grab_inds)
					# print (i,j)
					all_lons = lon[grab_inds]
					all_lats = lat[grab_inds]
					all_hli = h_li[grab_inds]
					# mean_lat = np.mean(all_lats)
					

					all_b_conf = bsnow_conf[grab_inds]
					all_b_od = bsnow_od[grab_inds]
					all_b_h = bsnow_h[grab_inds]

					mean_time = np.mean(dt[grab_inds])
					mean_year = np.mean(dy[grab_inds])
					filt = np.where(all_b_h<10000)
					bconf_mean = np.mean(all_b_conf[filt])
					bod_mean = np.mean(all_b_od[filt])
					bh_mean = np.mean(all_b_h[filt])
					bconf_x2 = np.sum(all_b_conf[filt]**2)
					bod_x2 = np.sum(all_b_od[filt]**2)
					bh_x2 = np.sum(all_b_h[filt]**2)
					pc_seg = len(filt[0])/N*100
					
					hli_mean = np.mean(all_hli[filt])
					df2 = [lastr, all_lons.mean(), all_lats.mean(), bconf_mean, bconf_x2, bod_mean, bod_x2, bh_mean,
						   bh_x2, N, len(filt[0]), pc_seg, hli_mean, mean_time, mean_year]
					df.loc[l] = df2
					l+=1
				else:
					pass

	f.close()
	return df


def downsample_all(dir,output_filename):
	all_files = os.listdir(dir)
	if '.DS_Store' in all_files: all_files.remove('.DS_Store')
	remove_these = []
	for fn in all_files:
		if '.h5' not in fn:
			remove_these.append(fn)
	for rfn in remove_these:
		all_files.remove(rfn)
	
	all_df = []
	for file in all_files:
		all_df.append(downsample(dir + '/' + file))
	final_df = pd.concat(all_df)
	final_df.to_csv(dir + '/' + output_filename)
import h5py 
import numpy as np
import os, glob
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import matplotlib.path as mpath

def show_plot(ax1):
	ax1.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
	ax1.add_feature(cartopy.feature.LAND)
	ax1.add_feature(cartopy.feature.OCEAN)
	theta = np.linspace(0, 2*np.pi, 100)
	center, radius = [0.5, 0.5], 0.5
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)
	ax1.set_boundary(circle, transform=ax1.transAxes)
	# plt.show()

	return

def plot_tracks(lon, lat, variable=None, **kwargs):
	# global ax1
	if variable is not None:
		plt.scatter(lon, lat, c=variable, transform = ccrs.PlateCarree(), **kwargs)
	else:
		plt.plot(lon, lat, transform = ccrs.PlateCarree(), **kwargs)
	return 

def main():
	#################################################
	os.chdir("./Small_data/")
	ls = glob.glob("*.h5")
	f = h5py.File(ls[0], 'r') 
	#################################################

	laser_strs = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
	ax1 = plt.subplot(111, projection=ccrs.SouthPolarStereo())	

	for i in range(len(laser_strs)):
		lastr = laser_strs[i]
		lat = f[lastr+'/lat'][:]
		lon = f[lastr+'/lon'][:]
		bsnow_conf = f[lastr+'/bsnow_conf'][:]
		dt = f[lastr+'/t_sec'][:]
		h_li = f[lastr+'/h_li'][:]
		plot_tracks(lon, lat, bsnow_conf)
		plot_tracks(lon, lat, None)


	cbar = plt.colorbar()
	cbar.set_label("Colorbar Label")
	show_plot(ax1)
	plt.show()
	# ax1.pcolormesh(xx, yy, variables[j], transform=data_crs)

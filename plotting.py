import numpy as np
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import cartopy
from netCDF4 import Dataset
from datetime import datetime, timedelta
os.chdir('/home/jeffrey/Snowblower/')
# sys.path.append('/Users/Ethan/Documents/Research/Git/SnowBlower/')

dg = Dataset("/home/jeffrey/Snowblower/Data/ERA_Interim_processed/erai_SH_analysis.nc", "r")

lons = dg.variables['lons'][:]
lats = dg.variables['lats'][:]
time = dg.variables['time'][:]
sktemp  = dg.variables['skt'][:]
temp_2m = dg.variables['t2m'][:]
u10, v10 = dg.variables['u10'][:], dg.variables['v10'][:]
wind_speed = np.sqrt(u10**2+v10**2)

xx, yy = np.meshgrid(lons, lats)


cartopy.crs.SouthPolarStereo(central_longitude=0.0, true_scale_latitude=-80)
sk0 = sktemp[0]

variables = [sk0 , u10[0]]


fig = plt.figure(figsize=[10, 5])
for j in range(0, 2):
	ax1 = plt.subplot(1, 2, j+1, projection=ccrs.SouthPolarStereo())
	# fig.subplots_adjust(bottom=0.05, top=0.95,
	                    # left=0.04, right=0.95, wspace=0.02)

	ax1.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
	ax1.add_feature(cartopy.feature.LAND)
	ax1.add_feature(cartopy.feature.OCEAN)
	theta = np.linspace(0, 2*np.pi, 100)
	center, radius = [0.5, 0.5], 0.5
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)

	ax1.set_boundary(circle, transform=ax1.transAxes)
	data_crs = ccrs.PlateCarree()

	ax1.pcolormesh(xx, yy, variables[j], transform=data_crs)

plt.show()

dt = datetime(1900 ,1, 1)
time2 = [timedelta(hours=int(x)) + dt for x in time]


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


variables = [sktemp[0] , wind_speed[0]]

from bokeh.plotting import figure, show, output_file

# x = np.linspace(lons.min(), lons.max(), len(lons))
# y = np.linspace(lats.min(), lats.max(), len(lats))

# p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
# p.x_range.range_padding = p.y_range.range_padding = 0

# # must give a vector of image data for image parameter
# p.image(image=[sk0], x=-180, y=-90, dw=360, dh=40, palette="Viridis11")

# output_file("image.html", title=" example")

# show(p)  # open a browser

import xarray as xr, hvplot.xarray, cartopy.crs as crs, geoviews as gv



xr.open_dataset()
proj = ccrs.SouthPolarStereo()
data = [xx, yy, sktemp[0]]
gv.tile_sources.ESRI * data.hvplot.points('Longitude', 'Latitude', geo=True, color='red', alpha=0.2, height=500,  xlim=(-180, -30), ylim=(0, 72))
import pandas as pd
idx = pd.date_range('1/1/2000', periods=1000)
df  = pd.DataFrame(np.random.randn(1000, 4), index=idx, columns=list('ABCD')).cumsum()

import hvplot.pandas
df.hvplot()




fig = plt.figure(figsize=[10, 5])
for j in range(0, 2):
	ax1 = plt.subplot(1, 2, j+1, projection=ccrs.SouthPolarStereo())
	# fig.subplots_adjust(bottom=0.05, top=0.95,
	                    # left=0.04, right=0.95, wspace=0.02)

	ax1.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())
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
time2 = [timedelta(hours=int(i)) + dt for i in time]


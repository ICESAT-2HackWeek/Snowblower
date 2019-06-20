import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plot_tracks import show_plot, plot_tracks

fakelons = 180 - np.random.random(50000)*360
fakelats = (np.random.random(50000)*15+75)*-1
fake_var = (np.random.random(50000)*0.7+0.3)*fakelons



lon_range = np.arange(-180, 180, 0.75)
lat_range = np.arange(-60, -90, -0.75)

lon_inds = np.digitize(fakelons, lon_range)
lat_inds = np.digitize(fakelats, lat_range)

meaned_data = np.zeros((lon_range.size, lat_range.size))
std_data = np.zeros((lon_range.size, lat_range.size))
## indexed as [lon, lat]

for i in range(len(meaned_data)):
	for j in range(len(meaned_data[i])):
		
		grab_inds = np.argwhere((lat_inds == j)*(lon_inds==i))
		#do we want to set a test for nans? if indexes are empty
		meaned_data[i, j] = np.mean(fake_var[grab_inds])
		std_data[i, j] = np.std(fake_var[grab_inds])

fig = plt.figure(figsize=(11,8))
## below section is required to get a nice triangular arrangement
gs0 = fig.add_gridspec(2, 1)
gs00 = gs0[0].subgridspec(1, 6)
gs01 = gs0[1].subgridspec(1, 6)

ax1 = fig.add_subplot(gs00[0,1:5], projection=ccrs.SouthPolarStereo())
# .subplot(2,2,(1,2), projection=ccrs.SouthPolarStereo())	
plt.title("Raw data, with noise added")
plot_tracks(fakelons, fakelats, variable=fake_var, cmap='bwr')
cbar0 = plt.colorbar(ax=ax1)
cbar0.set_label("Raw longitude")
show_plot(ax1)

# ax2 = plt.subplot(223, projection=ccrs.SouthPolarStereo())	
ax2 = fig.add_subplot(gs01[0,:3], projection=ccrs.SouthPolarStereo())

plt.title("Binned (mean) data")
data_crs = ccrs.PlateCarree()
xx, yy = np.meshgrid(lon_range, lat_range)
plt.pcolormesh(xx, yy, meaned_data.T, transform=data_crs, cmap='bwr')
cbar = plt.colorbar(ax=ax2)
cbar.set_label("Colorbar Label (unit)")
show_plot(ax2)

# ax3 = plt.subplot(224, projection=ccrs.SouthPolarStereo())	
ax3 = fig.add_subplot(gs01[0,3:], projection=ccrs.SouthPolarStereo())
plt.pcolormesh(xx, yy, std_data.T, transform=data_crs)
cbar2 = plt.colorbar()
cbar2.set_label("Standard deviation")
show_plot(ax3)


plt.tight_layout()
plt.show()

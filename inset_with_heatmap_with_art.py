%matplotlib inline
from utilities.color_map import *
from utilities.map_features import *
from utilities.plotting import *
from utilities.polygon_selection import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.mlab import griddata
import matplotlib.colors as colors
from matplotlib.cbook import get_sample_data
from matplotlib.colors import LinearSegmentedColormap
from numba import jit

import copy

from utilities.polygon_selection import *

from numba import jit

from obspy.imaging.beachball import beach

import myCalc as calc
import myProjection as projection

from bs4 import BeautifulStoneSoup as Soup
slip_df = pd.read_csv('data.txt', delim_whitespace=True, names=['lat', 'lon', 'depth', 'slip'])
slip_df.columns = ['la', 'lo', 'depth', 'slip']

output = calc.projection(fm_main[0], fm_main[1], fm_main[2], slip_df.copy())
data = output.copy()
numcols, numrows = 30, 80
xi = np.linspace(data.xx.min(), data.xx.max(), numcols)
yi = np.linspace(data.depth.min(), data.depth.max(), numrows)
xi, yi = np.meshgrid(xi, yi)

#-- Interpolate at the points in xi, yi
# "griddata" expects "raw" numpy arrays, so we'll pass in
# data.x.values instead of just the pandas series data.x
x, y, z = data.xx.values, data.depth.values, data.slip.values
zi = griddata(x, y, z, xi, yi, interp='linear')

fig, ax1 = plt.subplots(1, figsize=(15,4))

ax1.set_xlim(-45, 25)
ax1.set_ylim(0, 25)

bounds = np.array(np.concatenate([np.arange(0.5, 1.0, 0.01),np.arange(1.0,1.6,0.028)]))
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

c_bval = ax1.pcolormesh(r, c, pt_masked.transpose(), vmin=0.5, vmax=1.6, cmap='jet_r', norm=norm)
c_slip = ax1.contour(xi, yi, zi, colors='black', levels=range(1,8), alpha=0.75)
c_slip = ax1.contour(xi, yi, zi, cmap=cmap2, levels=range(1,8), alpha=0.95)
ax1.clabel(c_slip, inline=1, fontsize=10)

ticks = list(np.arange(0.5, 1.0, 0.1))
ticks.append(1)
ticks.append(1.5)
fig.colorbar(c_bval, ticks=ticks, label='b value', pad=0.01)

ax1.plot([0, 0], [0, 25], linestyle='--', linewidth=2, color='black', alpha=0.5)
ax1.plot(0, 12, marker='*', markersize=20, color='yellow')
ax1.plot([-23.455671, -23.455671], [0, 25], color='red', alpha=0.5, linestyle='--', linewidth=2)

ax1.set_yticks(range(0,25,5))
# ax1.set_yticks(range(0, 25, 1))
# ax1.set_xticks(range(-60, 30, 2))
ax1.set_xticks(range(-45, 25, 5))


im = plt.imread(get_sample_data('/home/max/pymap/research/b-value vs slip and stress drop/volcano-clipart-mount-fuji.png'))

# fig, ax = plt.subplots()
# ax.plot(range(10))
newax = fig.add_axes([0.41, 0.87, 0.2,  0.2], anchor='NE', zorder=-1)
newax.imshow(im)
newax.axis('off')

ax1.set_xlabel('Along Strike (km)', fontsize=15)
ax1.set_ylabel('Depth (km)', fontsize=15)

ax1.annotate('Mw7.0 main shock', xy=(0.2,12.2), xytext=(15, 18)
       ,arrowprops=dict(facecolor='black', shrink=0.005
                       ,width=0.2, headwidth=8))

ax1.annotate('Aso Caldera', xy=(-23.455671, 0.1), xytext=(-35, 2.5)
           ,arrowprops=dict(facecolor='black', shrink=0.005
                           ,width=0.2, headwidth=8))


ax1.text(s='NE', x=-40, y=-0.51, fontsize=20)
ax1.text(s='SW', x=20, y=-0.51, fontsize=20)

# kmeq = kuma[kuma.xx.between(-60,60) & kuma.yy.between(-1,5) & (kuma.timestamp < '2016-04-15')].copy()
kmeq = data_for_fmd

# kmeq.plot(ax=ax1, kind='scatter', marker='.', x='xx', y='depth'
#             ,s=np.exp(kmeq.mag), alpha=0.25, color='gray')
data_for_fmd.plot(ax=ax1, kind='scatter', marker='.', x='xx', y='depth'
            ,s=np.exp(data_for_fmd.mag), alpha=0.25, color='gray')

# ax1.grid(True)
ax1.set_xlabel('Along Strike (km)', fontsize=15)
ax1.set_ylabel('Depth (km)', fontsize=15)

# ax2 = fig.add_axes([left, bottom, width, height])
########################
# INSET LOWER          #
########################
node = (-17, 12.5)
eqs = find_earthquakes_within_circle(earthquakes=data_for_fmd[['xx','depth','idx']].values
                                      , circle_origin=node
                                      , circle_radius=5)

ax2 = fig.add_axes([0.55, 0.00, 0.1, 0.25])
h, e = np.histogram(data_for_fmd.loc[eqs].mag.values, bins=110, range=(0,10))
ch = np.cumsum(h[::-1])
ax2.plot(e[:-1], h, marker='s', linestyle='', color='None')
ax2.plot(e[::-1][:-1], ch, marker='^', linestyle='None', color='red')

ax2.set_yscale('log')

ax2.set_xlabel('Magnitude')
ax2.set_xticks([0,5,10])
ax2.set_yticks([1e0,1e1,1e2])

from matplotlib.patches import ConnectionPatch

xyA = node
xyB = (10, 1e2)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA='data', coordsB='data', axesA=ax1, axesB=ax2, zorder=20)
ax1.add_artist(con)

xyB = (0, 1e0)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA='data', coordsB='data', axesA=ax1, axesB=ax2, zorder=20)
ax1.add_artist(con)

data = fmd_stats.copy()
data[:,0:5][data[:,3] <= 50] = np.nan
r, c, pt = pivot_table(rows=data[:,5], cols=data[:,6], data=data[:,1])
pt_masked = np.ma.masked_invalid(pt)
datatest = pd.DataFrame(data, columns=['a','b','bstd','n','mc','xx','y'])
datatest = datatest[(datatest.xx==xyA[0]) & (datatest.y==xyA[1])].copy()
s = 'n={n}\nb={b}\nmc={mc}'.format(n=int(datatest.n.values[0])
                                  ,b=round(datatest.b.values[0],2)
                                  ,mc=round(datatest.mc.values[0],2))
ax2.text(s=s, x=5, y=1e1)

########################
# INSET UPPER          #
########################
node = (-16, 5)
eqs = find_earthquakes_within_circle(earthquakes=data_for_fmd[['xx','depth','idx']].values
                                      , circle_origin=node
                                      , circle_radius=5)

ax2 = fig.add_axes([0.35, 0.95, 0.1, 0.25])
h, e = np.histogram(data_for_fmd.loc[eqs].mag.values, bins=110, range=(0,10))
ch = np.cumsum(h[::-1])
ax2.plot(e[:-1], h, marker='s', linestyle='', color='None')
ax2.plot(e[::-1][:-1], ch, marker='^', linestyle='None', color='red')

ax2.set_yscale('log')

# ax2.set_xlabel('Magnitude')
ax2.set_title('GR-distribution')
ax2.set_xticks([0,5,10])
ax2.set_yticks([1e0,1e1,1e2])

from matplotlib.patches import ConnectionPatch

xyA = node
xyB = (10, 1e2)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA='data', coordsB='data', axesA=ax1, axesB=ax2, zorder=20)
ax1.add_artist(con)

xyB = (0, 1e0)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA='data', coordsB='data', axesA=ax1, axesB=ax2, zorder=20)
ax1.add_artist(con)

data = fmd_stats.copy()
data[:,0:5][data[:,3] <= 50] = np.nan
r, c, pt = pivot_table(rows=data[:,5], cols=data[:,6], data=data[:,1])
pt_masked = np.ma.masked_invalid(pt)
datatest = pd.DataFrame(data, columns=['a','b','bstd','n','mc','xx','y'])
datatest = datatest[(datatest.xx==xyA[0]) & (datatest.y==xyA[1])].copy()
s = 'n={n}\nb={b}\nmc={mc}'.format(n=int(datatest.n.values[0])
                                  ,b=round(datatest.b.values[0],2)
                                  ,mc=round(datatest.mc.values[0],2))
ax2.text(s=s, x=5, y=1e1)

ax1.invert_yaxis()
ax1.invert_xaxis()

# fig.savefig('/home/me/figures/cool.pdf', bbox_inches='tight')

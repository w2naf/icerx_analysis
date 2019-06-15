#!/usr/bin/python3
import matplotlib as mpl
mpl.use('Agg')

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

font = {'weight' : 'bold', 'size'   : 22}
mpl.rc('font', **font)


#class LowerThresholdRobinson(ccrs.Robinson):
#    @property
#    def threshold(self):
#        return 1e3

locs        = {}
locs['CHU'] = l = {}
l['lat']    =  45.291165502 # CHU
l['lon']    = -75.753663652 # CHU
l['lbl']    = 'CHU'
l['ha']     = 'center'
l['va']     = 'bottom'

locs['McM'] = l = {}
l['lat']    = -77.80 # McMurdo Station
l['lon']    = 166.67 # McMurdo Station
l['lbl']    = 'Arrival Heights\nMcMurdo Station'
l['ha']     = 'center'
l['va']     = 'bottom'

locs['WO'] = l = {}
l['lat']    =  40.785628
l['lon']    = -74.226118
l['lbl']    = 'West Orange, NJ'
l['ha']     = 'left'
l['va']     = 'top'

projection  = ccrs.PlateCarree(central_longitude=-135.)

scl = 2
fig = plt.figure(figsize=(scl*10,scl*8))
ax  = fig.add_subplot(1,1,1, projection=projection)

ax.stock_img()



l_0 = locs['CHU']
l_1 = locs['McM']

ax.plot([l_0['lon'], l_1['lon']], [l_0['lat'], l_1['lat']],
         color='red', linewidth=4,transform=ccrs.Geodetic())

l_0 = locs['CHU']
l_1 = locs['WO']
ax.plot([l_0['lon'], l_1['lon']], [l_0['lat'], l_1['lat']],
         color='blue', linewidth=4,transform=ccrs.Geodetic())


for key,l in locs.items():
#    ax.text(l['lon']- 3, l['lat'] - 12, l['lbl'],
#             horizontalalignment='right',
#             transform=ccrs.Geodetic())

    ax.text(l['lon'], l['lat'], l['lbl'],
             ha=l.get('ha','right'),
             va=l.get('va','bottom'),
             transform=ccrs.Geodetic())

ax.coastlines()
xlocs = [-220,-180,-120,-60,-20]
ylocs = np.arange(-80,81,20)
ax.gridlines(draw_labels=True,xlocs=xlocs,ylocs=ylocs)

ax.set_extent([-220, -20, -90, 90], crs=ccrs.PlateCarree())
fig.savefig('map.png',bbox_inches='tight')
plt.close(fig)

#!/usr/bin/env python
#
# Plotting topography and erosion rate
#
# This script uses matplotlib under ipython's pylab environment

import sys, os
import numpy as np
from matplotlib import pyplot as plt
import flac

###############################################

frame = int(sys.argv[1])

fl = flac.Flac()
x, z = fl.read_mesh(frame)
nx, nz = fl.nx, fl.nz  # number of nodes

xx = x[:,0]  # in km
topo = z[:,0]  # in km


sizeoffloat = 4
offset = (frame-1) * sizeoffloat * nx
with open('dtopo.0') as f:
    f.seek(offset)
    erosion_rate = -np.fromfile(f, dtype=np.single, count=nx)  # in mm/yr

fig, ax1 = plt.subplots()
#fig.clf()

ax1.plot(xx, topo, 'b-')
ax1.set_xlabel('x (km)')
ax1.set_ylabel('topo (km)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(xx, erosion_rate, 'r--')
ax2.set_ylabel('erosion rate (mm/yr)', color='r')
ax2.tick_params('y', colors='r')

#fig.tight_layout()
plt.show()

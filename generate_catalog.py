"""
Calculate a synthetic catalog of events matching a desired Gutenberg-Richter 
relationship
"""

import gutenbergrichter as gr
import math
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as P
from tqdm import tqdm
import cPickle as pickle
import sys

# Determine if the output filename is specified on command line
if len(sys.argv) > 1:
    fileout = sys.argv[1]
else:
    fileout = 'catalog.pkl'

# Basic characteristics of seismicity catalog
m0total = 1.0e16
max_m0 = math.pow(10.0,18.0)
minM = -1.0
min_m0 = gr.calc_m0(minM)
slope = 1.0

# Define the Gutenberg-Richter relationship values
gr_obj = gr.GutenbergRichter(b=slope, m0total=m0total, max_m0=max_m0,
                             min_m0=min_m0)
gr_obj.calc_a()

# Make a plot of the chosen G-R relationship
maxM = gr.calc_Mw(max_m0)
Msamp = 0.25
Mws = np.arange(minM, maxM + Msamp, Msamp)
Ns = gr_obj.get_N(Mws)
max_dep = 2.0 #max depth in km, events will be uniform between 0 and max_dep

# Calculate probability of a given event per second
secday = 60.0*60.0*24.0
secyear = secday*365.0
secmonth = secday*30.0
Nsec = Ns/secyear

# Generate a couple months of catalog
catlength = 2.0*secday
#(catalog, Nsc, Mwsc) = gr_obj.generate_catalog(catlength)
#(catalog2, Nsc2, Mwsc2) = gr_obj.generate_catalog(secmonth)
gr_obj.generate_catalog(catlength, max_dep=max_dep)

# Plot it all up
try:
    time_id = gr_obj.catalog.id_dict['time']
    mag_id = gr_obj.catalog.id_dict['magnitude']
    delta_id = gr_obj.catalog.id_dict['delta']
except AttributeError:
    time_id = 0
    mag_id = 1
    delta_id = 2

plt.figure(figsize=(10,15))
plt.subplot(3, 1, 1)
plt.semilogy(Mws,Ns,gr_obj.catalog.Mws,gr_obj.catalog.Ns)
plt.title("Gutenberg-Richter relationship")
plt.xlabel("Mw")
plt.ylabel("N")

ndays = int(catlength/secday)
plt.subplot(3, 1, 2)
plt.scatter(gr_obj.catalog.data[:,time_id]/secday,
            gr_obj.catalog.data[:,mag_id])
plt.title("%d day catalog" % ndays)
plt.xlabel("Day")
plt.ylabel("Magnitude")
plt.xlim([0, ndays])

plt.subplot(3, 1, 3)
numBins = 30
plt.hist(gr_obj.catalog.data[:,delta_id],numBins,color='green')
plt.title("Distances")
plt.xlabel("Distance (degrees)")
plt.ylabel("Frequency")

#plt.show()
figname = 'catalog.png'
P.savefig(figname)

# Write out catalog to pickle file
filename = fileout
with open(filename, 'wb') as f:
    pickle.dump(gr_obj, f, -1)





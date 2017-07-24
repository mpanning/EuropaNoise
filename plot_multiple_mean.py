from obspy.signal import PPSD
from obspy import read, UTCDateTime
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
#from obspy.imaging import cm
import sys
import glob
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from operator import itemgetter


# Get all filenames with this wildcard
# Make a 4 panel figure with all 5 seismicity models for each Europa model
# comp = sys.argv[1]
comp = 'MXZ'
models = ['ice5', 'ice20', 'ice5_lowQ', 'ice20_lowQ'] # Later add low Q
catalogs = ['modelA', 'modelB', 'pref', 'modelC', 'modelD']
labels = ['Model A', 'Model B', 'Preferred', 'Model C', 'Model D']
# catalogs = ['modelA', 'modelB', 'modelC', 'modelD']
# labels = ['Model A', 'Model B', 'Model C', 'Model D']
linecolors = ['b', 'g', 'm', 'r', 'c'] 
titles = ['5 km ice shell', '20 km ice shell', '5 km ice shell, low Q',
          '20 km ice shell, low Q']


# Define a flat response file and other variables needed
paz = {'gain': 1.0,
       'poles': [],
       'zeros': [],
       'sensitivity': 1.0}
secday = 24.0 * 3600.0

# Load SP self-noise data
seis_types = ['SP_Imperial', 'STS2', '10Hz_geophone',
              'Trillium_compact']

period = dict()
psd = dict()
style = dict()
size = dict()
color = dict()
seislabels = dict()
for seis in seis_types:
    filename = 'noise_%s.txt' % seis
    print('Loading %s' % filename)
    selfnoise = np.loadtxt(filename)
    period[seis] = 1. / (selfnoise[1:, 0])
    psd[seis] = 10 * np.log10(selfnoise[1:, 1]) * 2
    if seis == 'SP_Imperial':
        color[seis] = 'darkgreen'
        style[seis] = '--'
        size[seis] = 2.0
    else:
        style[seis] = '--'
        size[seis] = 2.0

color['STS2'] = 'darkred'
color['10Hz_geophone'] = 'black'
color['Trillium_compact'] = 'fuchsia'
seislabels['SP_Imperial'] = 'SP'
seislabels['STS2'] = 'STS2'
seislabels['10Hz_geophone'] = 'geophone'
seislabels['Trillium_compact'] = 'T-compact'

# Loop over models to make figure panels
fig = plt.figure(figsize=(13.0, 13.0))
for i, model in enumerate(models):
    # Loop over catalogs to make lines in figure panel
    for j, catalog in enumerate(catalogs):
        db_short = model + '.' + catalog + '_cat'
        filenames = glob.glob('noise_records/*%s*.%s' % (db_short,comp))

        print 'noise_records/*%s*.%s' % (db_short,comp)
        filenames.sort()

        print filenames

        st = read(filenames[0])
        # Differentiate to velocity, which is what PPSD expects
        st.differentiate()
        tr_tmp = st[0]
        tr_tmp.stats.channel = comp

        ppsd = PPSD(tr_tmp.stats, paz, db_bins=[-300, -50, 5],
                    period_limits=[0.5, 500])
        ndays = int(math.ceil((tr_tmp.stats.endtime
                               - tr_tmp.stats.starttime)/secday))

        # Load all files in list filenames and modify the startdates to make
        # them consecutive
        iday = 0
        for filename in filenames:
            iday += ndays
            st = read(filename)
            # Differentiate to velocity
            st.differentiate()
            st[0].stats.starttime = UTCDateTime(2028, 1, iday, iday, 0, 0)
            st[0].stats.channel = comp
            ppsd.add(st)


            # Plot PPSD
            # h = ppsd.plot(period_lim=[0.5,500], show=False,
            #              show_coverage=False, max_percentage=15)
            # plt.figure(num=1, figsize=(20.0,20.0))
            # ppsd.plot(period_lim=[0.5,500], show=False,
            #           show_coverage=False, max_percentage=12)
            # ax = h.axes[0]
        ax = plt.subplot(2, 2, i+1)

        (meanpd, meanpsd) = ppsd.get_mean()
        ax.semilogx(meanpd, meanpsd, linecolors[j], linewidth=3,
                    label=labels[j])
        plt.xlim([0.5, 500])
        plt.ylim([-300, -50])

    # Add SP self noise plot
    for seis in seis_types:
        line2, = ax.plot(period[seis], psd[seis], linewidth=size[seis],
                        color=color[seis], linestyle=style[seis],
                        label=seislabels[seis])
    
    # Plot Earth noise models
    nlnm_pd, nlnm = get_nlnm()
    nhnm_pd, nhnm = get_nhnm()
    ax.plot(nhnm_pd, nhnm, linewidth=2, color='darkgrey',
            label='Earth noise')
    ax.plot(nlnm_pd, nlnm, linewidth=2, color='darkgrey')

    # Add axis labels
    if i > 1:
        plt.xlabel('Period (s)')
    if (i == 0 or i == 2):
        plt.ylabel(r'Power ($m^2/s^4/Hz$) (dB)') 

    plt.title(titles[i])
    if i == 1:
        mpl.rcParams['legend.fontsize'] = 'small'
        ax.legend(bbox_to_anchor=(1.1, 1), fancybox=True)


# Add in panel labels
plt.annotate('A', xy=(0.15, 0.87), xycoords='figure fraction', size='x-large')
plt.annotate('B', xy=(0.57, 0.87), xycoords='figure fraction', size='x-large')
plt.annotate('C', xy=(0.15, 0.43), xycoords='figure fraction', size='x-large')
plt.annotate('D', xy=(0.57, 0.43), xycoords='figure fraction', size='x-large')
    
plt.savefig('Summary_noise.png')

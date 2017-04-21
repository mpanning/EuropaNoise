from obspy.signal import PPSD
#from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
from obspy import read, UTCDateTime
#from obspy.imaging import cm
import sys
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from operator import itemgetter


# Get all filenames with this wildcard
comp = sys.argv[1]
db_short = sys.argv[2]
filenames = glob.glob('noise_records/*%s*.%s' % (db_short,comp))

print 'noise_records/*%s*.%s' % (db_short,comp)
filenames.sort()

print filenames

# Define a response function that translates the
# displacement seismograms into velocity (which is what
# PPSD uses internally.
# paz = {'gain': 1.0,
#        'poles': [],
#        'zeros': [0j],
#        'sensitivity': 1.0}

# Define a flat response file
paz = {'gain': 1.0,
       'poles': [],
       'zeros': [],
       'sensitivity': 1.0}

# tr_tmp = read(filenames[0])[0]
st = read(filenames[0])
# Differentiate to velocity, which is what PPSD expects
st.differentiate()
tr_tmp = st[0]
tr_tmp.stats.channel = comp

# Calculate psd using mlab.psd from matplotlib
dt = st[0].stats['delta']
fs = 1. / dt
length = len(st[0].data)
# differentiate one more time to acceleration
st.differentiate()
# Get the psd, but throw out the plot
# fig = plt.figure()
# (Pxx, freqs) = plt.psd(st[0].data, NFFT=length//16, pad_to=length, Fs=fs,
#                       sides='onesided', scale_by_freq=True)
# plt.close(fig) 
# Define nfft and overlap identical to obspy.signal.PPSD
ppsd_length = 3600.0
nfft = ppsd_length * fs
nfft = nfft / 4.0
nfft = int(math.pow(2, math.floor(math.log(nfft, 2))))

# Try a range of nfft values for mlab.psd
# nffts = np.array([nfft, ppsd_length * fs, length])
nffts = np.array([nfft])
nffts = nffts.astype(int)
nlaps = 0.75 * nffts
nlaps = nlaps.astype(int)
psd_color = 'lime'
psd_label = 'Highest amplitude segments'
npsd = 3
Ps = []
freqs = []
power = []
for i in range(len(nffts)):
    # (Pxx, freq) = mlab.psd(st[0].data, NFFT=nffts[i], Fs=fs,
    #                        sides='onesided', scale_by_freq=True,
    #                        detrend=mlab.detrend_linear, noverlap=nlaps[i])
    nslices = int(math.floor(length/ppsd_length))
    for j in range(nslices):
        dataslice = st[0].data[j*nffts[i]:(j+1)*nffts[i]]
        (Pxx, freq) = mlab.psd(dataslice, NFFT=nffts[i], Fs=fs,
                               sides='onesided', scale_by_freq=True,
                               detrend=mlab.detrend_linear, noverlap=nlaps[i])
    
        # Remove DC term and convert to db
        Pxx = 10.0 * np.log10(Pxx[1:])
        power.append(Pxx.sum())
        freq = freq[1:]

        Ps.append(Pxx)
        freqs.append(freq)

# Sort by summed power
zipsort = sorted(zip(power, Ps, freqs), key=itemgetter(0))
Ps = [x[1] for x in zipsort]
freqs = [x[2] for x in zipsort] 

ppsd = PPSD(tr_tmp.stats, paz, db_bins=[-300, -50, 5],
            period_limits=[0.5, 500])
secday = 24.0 * 3600.0
ndays = math.ceil((tr_tmp.stats.endtime - tr_tmp.stats.starttime)/secday)

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

# Load SP self-noise data
seis_types = ['SP_Imperial', 'STS2', '10Hz_geophone',
              'Trillium_compact']
period = dict()
psd = dict()
style = dict()
size = dict()
color = dict()
for seis in seis_types:
    filename = 'noise_%s.txt' % seis
    print('Loading %s' % filename)
    selfnoise = np.loadtxt(filename)
    period[seis] = 1. / (selfnoise[1:, 0])
    psd[seis] = 10 * np.log10(selfnoise[1:, 1]) * 2
    if seis == 'SP_Imperial':
        color[seis] = 'r'
        style[seis] = '-'
        size[seis] = 2.5
    else:
        style[seis] = '--'
        size[seis] = 2.0

color['STS2'] = 'cyan'
color['10Hz_geophone'] = 'yellow'
color['Trillium_compact'] = 'fuchsia'

# Plot PPSD
#h = ppsd.plot(period_lim=[0.5,500], show=False,
#              show_coverage=False, max_percentage=15)
# plt.figure(num=1, figsize=(20.0,20.0))
ppsd.plot(period_lim=[0.5,500], show=False,
          show_coverage=False, max_percentage=12)
#ax = h.axes[0]
ax = plt.gca()
fig = plt.gcf()
fig.set_size_inches(10.0, 10.0, forward=True)

# Add in psd lines
# Pxx, freqs calculated above
# Plot npsd largest segments
Pdb = Ps[-npsd:]
freq = freqs[-npsd:]
for i in range(len(Pdb)):
    psd_pd = 1.0/freq[i]
    ax.plot(psd_pd, Pdb[i], color=psd_color, linestyle = '-', linewidth=0.5)
ax.plot(0., 0., color=psd_color, label=psd_label, linewidth=0.5)

# Add low and high Earth noise model
# nlnm_pd, nlnm = get_nlnm()
# nhnm_pd, nhnm = get_nhnm()
# ax.plot(nhnm_pd, nhnm, linewidth=2, color='darkgrey',
#         label='Earth noise model')
# ax.plot(nlnm_pd, nlnm, linewidth=2, color='darkgrey')

(meanpd, meanpsd) = ppsd.get_mean()
ax.plot(meanpd, meanpsd, linewidth=2, color='black', label='Mean PSD')
ax.plot(0., 0., linewidth=2, color='darkgrey', label='Earth noise model')
ax.legend(loc=1)

plt.savefig('noisemodel_%s_%s.png' % (db_short,comp))

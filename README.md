# EuropaNoise

All noise records generated for the paper "The seismic noise environment of Europa" are included in the subdirectory noise_records.

The noise records are all SAC files (Seismic Analysis Code http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/). 

The naming scheme of the files all are set up as `model`.`catalog`.`component`

The possible values for `model` are:
* ice5 - 5 km thick ice shell with high Q
* ice5_lowQ - ice5 with ice shell Q reduced by a factor of 10
* ice20 - 20 km thick ice shell with high Q
* ice20_lowQ - ice20 with ice shell Q reduced by a factor of 10
* ice5_scat - ice5 with heterogeneities as described in paper

The catalogs are named consistent with the seismicity model (models A through D or pref for preferred), as well as the specific catalog file as stored as pickle files in the catalogs subdirectory

`component` is either MXZ, MXE, or MXN for the vertical, east and north components of motion respectively.

The catalogs subdirectory contains pickle files with event catalogs stored as gutenbergrichter Catalog objects as defined in `gutenbergrichter.py`

Sample python scripts to generate catalogs, noise records and some figures from the paper are also included.  Currently the waveform databases are available upon request from the authors, but they are planned to be archived in a web accessible format in the near future.  Updates coming.

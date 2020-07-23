# MWISPdbscan
use DBSCAN to find molecular clouds and extract catalog

# NOTE:
Most of the files are useless and they are copied here only to make the script runable.

Author : Qing-Zeng Yan

Paper Link:  https://ui.adsabs.harvard.edu/abs/2020arXiv200613654Y



# V1.
example.py

Please make sure the example.py file run correctly

You may need to install lots of python packages.



Cloud Table colnames:
"_idx": index of clouds
"area_exact": Area of cloud
"v_cen": average radial veloicty, km/s
"v_rms":  weighted (by temperature) standard deviation of velocity, km/s
"x_cen": average Galactic longitude in degree
"l_rms":  same as "v_rms" but for the Galadctic longtitude in degree
"y_cen": average Galactic latitude in degree
"b_rms":  same as "v_rms" but for the Galadctic latitude in degree
"sum": sum of voxel values in Kelvin
"pixN": total number of voxels
"peak": Peak value of clouds
"peakL": Index of peak in PPV
"peakB": Index of peak in PPV
"peakV": Index of peak in PPV
"area_accurate": Exact area of cloud (considering the factor of cos(b), current not calculated for the speed purpose)
"allChannel": total number of channels
"hass22": has a 2x2 compact region (a beam)
 

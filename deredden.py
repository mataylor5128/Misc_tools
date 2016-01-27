def get_red(ra,dec,filter):

	do_ra = ra#201.365062792
	do_dec = dec#-43.0191125833

	# print cena_ra, cena_dec

	ra1 = int(do_ra/15.)
	ra2 = int((do_ra/15.-ra1)*60.)
	ra3 = ((do_ra/15.-ra1)*60.-ra2)*60.

	dec1 = int(dec)
	dec2 = np.abs(int((do_dec-dec1)*60.))
	dec3 = (np.abs((do_dec-dec1)*60.)-dec2)*60.

	do_filter = filter

	tempfile = 'temp.txt'
	f = open(tempfile,'w')
	call(['curl',
		'http://ned.ipac.caltech.edu/cgi-bin/calc?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000&lon=%i+%i+%f&lat=%i+%i+%f&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0' % (ra1,ra2,ra3,dec1,dec2,dec3)],
		stdout=f
		)
	f.close()
	f = open(tempfile,'r')
	data = f.read()
	red_value = float(data[data.find('SDSS    %s' % do_filter)+18:data.find('SDSS    %s' % do_filter,)+23])
	f.close()
	os.remove(tempfile)

	return red_value

from astropy.io import fits	
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
import os

do_filt = 'z'

hdulist = fits.open('survey_tile1_%s_psf.ldac' % do_filt)

do_filt = 'g'

tbldata = hdulist[2].data

ra_vals = []
dec_vals = []
red_vals = []
for ii in range(500):
	ind = np.random.random_integers(0,len(tbldata),1)
	temp_ra = tbldata['ALPHA_J2000'][ind][0]
	temp_dec = tbldata['DELTA_J2000'][ind][0]
	ra_vals.append(temp_ra)
	dec_vals.append(temp_dec)
	red_vals.append(get_red(temp_ra,temp_dec,do_filt))

cm = plt.cm.get_cmap('jet')

fig = plt.figure()
red_fig = plt.scatter(ra_vals,dec_vals,marker='o',edgecolor="None",
	c=red_vals,
	vmin=np.min(red_vals),
	vmax=np.max(red_vals),
	cmap=cm,s=5)

plt.xlim(np.max(tbldata['ALPHA_J2000']),np.min(tbldata['ALPHA_J2000']))
plt.ylim(np.min(tbldata['DELTA_J2000']),np.max(tbldata['DELTA_J2000']))
plt.xlabel(r'$\delta$ (J2000)',fontsize=18)
plt.ylabel(r'$\alpha$ (J2000)',fontsize=18)

cb = fig.colorbar(red_fig, orientation='vertical',pad = 0.05)#, aspect=20, orientation='vertical',pad = 0.05)
cb.set_label(r"%s'-band Extinction [mag]" % do_filt,fontsize=18)

plt.savefig('%s-band_red_map.pdf' % do_filt)
plt.show()
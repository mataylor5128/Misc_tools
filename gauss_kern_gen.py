def fwhm_convert(xx):
	return xx/(2*np.sqrt(2*np.log(2)))

import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian2DKernel
import sys

fwhm = sys.argv[1]
dim = sys.argv[2]
filter = sys.argv[3]
tile = sys.argv[4]

kernel = Gaussian2DKernel(fwhm_convert(float(fwhm)),x_size=int(dim),y_size=int(dim),mode='integrate')

# for line in kernel.array:
# 	print str(line)[1:-1].replace('\n','')
# exit()

kernel = np.sum(kernel.array)*kernel.array/np.max(kernel.array)

#plt.figure()
#plt.imshow(kernel,interpolation='none',origin='lower')
#plt.colorbar()
#plt.show()

#print 'gauss_'+fwhm+'_'+dim+'x'+dim+'_t'+tile+'_'+filter+'.conv'
fileout = open('gauss_'+fwhm+'_'+dim+'x'+dim+'_t'+tile+'_'+filter+'.conv','w')
print >> fileout, 'CONV NORM'
print >> fileout, '# '+dim+'x'+dim+' convolution mask of a gaussian PSF with FWHM = '+fwhm+' pixels.'
for line in kernel:
#	print str(line)[1:-1]
	print >> fileout, str(line)[1:-1].replace('\n','')

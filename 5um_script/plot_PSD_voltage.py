import numpy, h5py, matplotlib
import matplotlib.pyplot as plt
import os
import scipy.signal as sp
import numpy as np
import bead_util as bu
import glob

path = r"F:\data\20220317\5um_SiO2\3\test\plot_xy_pin\2"
fname0 = r"LPmbar_xyzcool_0.h5"
fname1 = r"LPmbar_xyzcool_1.h5"
fname2 = r"LPmbar_xyzcool_2.h5"
fname3 = r"LPmbar_xyzcool_3.h5"
realcsdnorm = False


Fs = 10e3  #
# this is ignored with HDF5 files
NFFT = 2**16

def getdata(fname):
    print ("Opening file: ", fname)
    ## guess at file type from extension
    _, fext = os.path.splitext( fname )
    if( fext == ".h5"):
        f = h5py.File(fname,'r')
        dset = f['beads/data/pos_data']
        dat = numpy.transpose(dset)
		#max_volt = dset.attrs['max_volt']
		#nbit = dset.attrs['nbit']
        Fs = dset.attrs['Fsamp']
        pid = dset.attrs['PID']		
		#dat = 1.0*dat*max_volt/nbit
        dat = dat * 10./(2**15 - 1)
        Press = dset.attrs['pressures'][0]
        print (Press)
    else:
	    dat = numpy.loadtxt(fname, skiprows = 5, usecols = [2, 3, 4, 5, 6])
    xpsd, freqs = matplotlib.mlab.psd(dat[:, bu.xi]-numpy.mean(dat[:, bu.xi]), Fs = Fs, NFFT = NFFT) 
    ypsd, freqs = matplotlib.mlab.psd(dat[:, bu.yi]-numpy.mean(dat[:, bu.yi]), Fs = Fs, NFFT = NFFT)
    zpsd, freqs = matplotlib.mlab.psd(dat[:, bu.zi]-numpy.mean(dat[:, bu.zi]), Fs = Fs, NFFT = NFFT)

    xpsd_outloop, freqs = matplotlib.mlab.psd(dat[:, 4]-numpy.mean(dat[:, 4]), Fs = Fs, NFFT = NFFT)
    # Ddrive = dat[:, bu.drive]*np.gradient(dat[:,bu.drive])
    # DdrivePSD, freqs =  matplotlib.mlab.psd(Ddrive-numpy.mean(Ddrive), Fs = Fs, NFFT = NFFT))
    print (pid)

    # f, Pxy = sp.csd(dat[:, 0]-numpy.mean(dat[:, 0]), dat[:, 4] - numpy.mean(dat[:, 4]), Fs, nperseg=NFFT)
    # plt.figure()
    # plt.loglog(f, np.sqrt(np.abs(np.real(Pxy))))
	# norm = numpy.median(dat[:, bu.zi])
    if realcsdnorm:
        f, Pxy = np.abs(np.real(sp.csd(dat[:, 0]-numpy.mean(dat[:, 0]), dat[:, 4] - numpy.mean(dat[:, 4]), Fs, nperseg=NFFT, scaling = "spectrum")))
        f, Pxx = sp.csd(dat[:, 0]-numpy.mean(dat[:, 0]), dat[:, 0]-numpy.mean(dat[:, 0]), Fs, nperseg=NFFT, scaling = "spectrum")
        f, Pyy = sp.csd(dat[:, 5]-numpy.mean(dat[:, 5]), dat[:, 5]-numpy.mean(dat[:, 5]), Fs, nperseg=NFFT, scaling = "spectrum")
        Cxy = (Pxy**2)/(Pxx*Pyy)
        return [freqs, xpsd, ypsd, dat, zpsd, xpsd_outloop, f, Cxy]
    return [freqs, xpsd, ypsd, dat, zpsd, xpsd_outloop]


data0 = getdata(os.path.join(path, fname0))
data1 = getdata(os.path.join(path, fname1))
data2 = getdata(os.path.join(path, fname2))
data3 = getdata(os.path.join(path, fname3))
a=np.where(data0[0]>10)[0][0]
b=np.where(data0[0]<50)[0][-1]
print(a,b)
Fs = 10000
fig = plt.figure()
plt.subplot(2, 2, 1)
plt.loglog(data0[0][a:b],  np.sqrt(data0[1][a:b]),label="no_drive_x")
plt.loglog(data0[0][a:b], np.sqrt(data0[2][a:b]),label="no_drive_y")
plt.legend()
plt.ylabel("V/rtHz")
plt.subplot(2, 2, 2)
plt.loglog(data1[0][a:b], np.sqrt(data1[1][a:b]),label="x_drive_x")
plt.loglog(data1[0][a:b], np.sqrt(data1[2][a:b]),label="x_drive_y")
plt.legend()
plt.subplot(2, 2, 3)
plt.loglog(data2[0][a:b], np.sqrt(data2[1][a:b]),label="y_drive_x")
plt.loglog(data2[0][a:b], np.sqrt(data2[2][a:b]),label="y_drive_y")
plt.legend()
plt.ylabel("V/rtHz")
plt.xlabel("Frequency[Hz]")
plt.subplot(2, 2, 4)
plt.loglog(data3[0][a:b],  np.sqrt(data3[1][a:b]),label="no_drive_x")
plt.loglog(data3[0][a:b], np.sqrt(data3[2][a:b]),label="no_drive_y" )
plt.xlabel("Frequency[Hz]")
plt.legend()
plt.show()


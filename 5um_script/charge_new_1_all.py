
## load all files in a directory and plot the correlation of the resonse
## with the drive signal versus time

import numpy as np
import matplotlib, calendar
import matplotlib.pyplot as plt
import os, re, time, glob
import bead_util as bu
import scipy.signal as sp
import scipy.optimize as opt
from scipy.fft import fft, ifft, rfft, irfft, fftfreq
import cPickle as pickle
import time

path = r"F:\data\20220302\5um_SiO2\3\discharge\1"
ts = 1.

fdrive = 35. #31.
make_plot = True

data_columns = [0, bu.xi] # column to calculate the correlation against
drive_column = 3 # column containing drive signal

def getdata(fname):
    print ("Processing ", fname)
    dat, attribs, cf = bu.getdata(os.path.join(path, fname))

    if( len(attribs) > 0 ):
        fsamp = attribs["Fsamp"]
    print ("Getting data from: ", fname) 
    dat, attribs, cf = bu.getdata(os.path.join(path, fname))
    fsamp = attribs["Fsamp"]
    xdat = dat[:,data_columns[1]]# data of the displacement
    xdat = xdat - np.mean(xdat)
    fourier_x=rfft(xdat)#Do the Fourier transform of the x
    f = fftfreq(len(xdat),1/fsamp)
    freq_x = np.abs(f[0:int(len(xdat)/2)+1])
        
    data_drive=dat[:,drive_column] - np.mean(dat[:,drive_column]) # Electric field without offset
    fourier_drive=rfft(data_drive) # Do the Fourier transform of the drive signal
    f = fftfreq(len(data_drive),1/fsamp)
    freq_drive = np.abs(f[0:int(len(data_drive)/2)+1])#frequency domain axis
    E_freq=freq_drive[np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0]]
    #plt.plot(freq_drive,fourier_drive) ## check if the Fourier transform is done as assumed
    #plt.show()
    print("Electric field's frequency is:", E_freq[0])
    #print(len(xdat),len(data_drive))
    Conv=fourier_drive*np.conjugate(fourier_x) 
    Corr=irfft(Conv)
    print(len(Corr))
    n=np.where(np.abs(freq_x-E_freq)<2)
    print(n,np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0][0])
    maxv = max(np.real(Conv[n]))
    print("Displacement peak is:",maxv)
    print("drive peak is:", max(np.abs(fourier_drive)))
    Phase=np.angle(fourier_x[np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0][0]])-np.angle(fourier_drive[np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0]])
    print("the phase is:", Phase[0])
    cf.close()
    return Phase[0], maxv

def get_most_recent_file(p):

    ## only consider single frequency files, not chirps
    filelist = glob.glob(os.path.join(p,"*.h5"))  ##os.listdir(p)
    #filelist = [filelist[0]]
    mtime = 0
    mrf = ""
    for fin in filelist:
        if( fin[-3:] != ".h5" ):
            continue
        f = os.path.join(path, fin) 
        if os.path.getmtime(f)>mtime:
            mrf = f
            mtime = os.path.getmtime(f)

    fnum = re.findall('\d+.h5', mrf)[0][:-3]
    return mrf#.replace(fnum, str(int(fnum)-1))


corr_data=[]
if make_plot:
    fig0 = plt.figure()
    # plt.hold(False)

last_file = ""
filelist = glob.glob(os.path.join(path, "*.h5"))
filelist = sorted(filelist, key=os.path.getmtime)
n=0

for cfile in filelist:
    n=n+1
    print(n,'s')
    corr = getdata(cfile)
    corr_data.append(corr)

if make_plot:
    plt.clf()
    plt.plot(np.array(corr_data))
    plt.grid()
    plt.xlabel('t(s)')
    plt.ylabel('Correlation value')
    plt.show()
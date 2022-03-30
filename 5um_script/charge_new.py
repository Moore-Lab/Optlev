
## load all files in a directory and plot the correlation of the resonse
## with the drive signal versus time

import numpy as np
import matplotlib, calendar
import matplotlib.pyplot as plt
import os, re, time, glob
import bead_util as bu
import scipy.signal as sp
import scipy.optimize as opt
from scipy.fft import fft, ifft, rfft, fftfreq
import cPickle as pickle
import time

path = r"F:\data\20220318\5um_SiO2\2\discharge\2"
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
    n=np.where(np.abs(freq_x-E_freq)<2)
    print(n,np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0][0])
    maxv = max(np.real(np.abs(fourier_x[n])))
    print("Displacement peak is:",maxv)
    print("drive peak is:", max(np.abs(fourier_drive)))
    Phase=np.angle(fourier_x[np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0][0]])-np.angle(fourier_drive[np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0]])
    print("the phase is:", Phase[0])
    #maxf=np.where(np.abs(fourier_x)==maxv)[0]
    #print(maxf)
    #print(fourier_x[np.where(np.abs(fourier_drive)==max(np.abs(fourier_drive)))[0]][0], np.real(maxv))
    cf.close()
    return Phase[0], np.real(maxv)

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
while( True ):
    ## get the most recent file in the directory and calculate the correlation

    cfile = get_most_recent_file( path )
    
    ## wait a sufficient amount of time to ensure the file is closed
    print (cfile)
    time.sleep(ts)

    if( cfile == last_file ): 
        continue
    else:
        last_file = cfile

    ## this ensures that the file is closed before we try to read it
    time.sleep( 1 )

    start_1 = time.time()
    corr = getdata( cfile )
    corr_data.append(corr )
    end_1 = time.time()
    print("time1=",end_1 - start_1)

    np.savetxt( os.path.join(path, "current_corr.txt"), [corr,] )
    if make_plot:
        plt.clf()
        start_2 = time.time()
        plt.plot(np.array(corr_data),label=["Phase","Charge"])
        plt.grid()
        plt.legend()
        plt.draw()
        end_2 = time.time()
        print("time2=",end_2 - start_2)
        plt.pause(0.01)
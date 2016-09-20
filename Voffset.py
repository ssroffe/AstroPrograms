import tools as tls
from astropy.io import fits
import numpy as np
import scipy.interpolate as sp
import matplotlib.pyplot as plt
from plot_format import plot_format
import os
import sys

lines = [line.rstrip('\n') for line in open('filelist')]
plot_format()
#print lines
plotNumber = 1
plotIndex = 1

for j in range(len(lines)):
    #plotIndex += 1
    path = lines[j]

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    timeTaken = basename[15:]
    #timeTaken = basename[7:]
    #print timeTaken

    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    
    Owl,Normflux,errorNorm = tls.NormNoPlot(path)
    
    ## Amplify the lines
    #Normflux = ((Normflux-1)*5)**3+1

    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10
    lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    lineWindows = np.array([[4060.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"])

    for lineIndex in range(1,4):
        
        #offset = 30
        offset = 15

        upperLine = lineList[lineIndex] + offset
        lowerLine = lineList[lineIndex] - offset
        
        #plt.axvline(upperLine,color='black')
        #plt.axvline(lowerLine,color="black")
        #plt.plot(Owl,Normflux)
        #plt.show()
        
        wherr = np.where((Owl >= lowerLine) & (Owl <= upperLine))
        flux = Normflux[wherr]
        wl = np.linspace(lowerLine,upperLine,len(flux))
        
        #plt.plot(wl,flux)
        #plt.show()
            
        vel = []
        for w in range(len(wl)):
            v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
            vel.append(v)
                
                
        #Plot Hbeta vertical line
        plt.subplot(2,2,plotIndex)
        plt.plot(Owl,Normflux)
        plt.axhline(1,ls='--',color='black')
        plt.subplot(2,2,plotIndex+1)
        plt.axvline(0,ls='--',color='black')
        
        plt.step(vel,flux+lineIndex*20.,where='mid',linewidth=2,label=lineNames[lineIndex])
        plt.xlim(min(vel),max(vel))
        plt.title("Spectrum: "+str(j))
        plt.xlabel("velocity [km/s]")
        plt.ylabel("Normalized flux + offset")
        plt.legend(prop={'size':5})

    plotIndex += 2    
    if plotIndex == 5:
        plt.savefig("vPlots/"+wdName+"_pg"+str(plotNumber)+".pdf")
        plot_format()
        plotNumber += 1
        plotIndex = 1


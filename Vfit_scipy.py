import tools as tls
from astropy.io import fits
import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt
from plot_format import plot_format
import os
import sys
import emcee as mc


def lorentzian(x, depth, width, RVShift):
    return 1.0-depth/(1.0 + ((x-RVShift)/width)**2)

lines = [line.rstrip('\n') for line in open('filelist')]
plot_format()
#print lines
plotNumber = 1
plotIndex = 1

for j in range(len(lines)):

    rvs = None

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
    #Normflux = ((Normflux-1)*2)+1

    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10
    lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    lineWindows = np.array([[4060.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"])
    #lineIndex = 1
    for lineIndex in range(1,4):
        
        #offset = 30
        offset = 25

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

        ############################
        ###### FITTING STUFF #######
        ############################
        p = [7.0,50,0.0]

        if rvs == None:
            popt, pcov = sp.curve_fit(lorentzian,vel,flux,maxfev=1500)
            rvs = popt[2]
        else:
            popt, pcov = sp.curve_fit(lambda x,depth,width: lorentzian(x,depth,width,rvs),vel,flux,maxfev=1500)

        #Plot Hbeta vertical line
        plt.subplot(2,2,plotIndex)
        plt.plot(Owl,Normflux)
        plt.axhline(1,ls='--',color='black')
        plt.subplot(2,2,plotIndex+1)
        plt.axvline(0,ls='--',color='black')
        plt.axvline(rvs,ls='--',color='purple',linewidth=1.5)
        off = 0.3

        plt.step(vel,flux+lineIndex*off,where='mid',linewidth=2,label=lineNames[lineIndex])
        plt.plot(vel,lorentzian(vel,popt[0],popt[1],rvs)+lineIndex*off,color='black',linewidth=2)
        plt.xlim(min(vel),max(vel))
        plt.title("Spectrum: "+str(j))
        plt.xlabel("velocity [km/s]")
        plt.ylabel("Normalized flux + offset")
        plt.legend(prop={'size':5})
        #plt.show()
    
    plotIndex += 2    
    if plotIndex == 5:
        plt.savefig("fitPlots/"+wdName+"_fit_pg"+str(plotNumber)+".pdf")
        plot_format()
        plotNumber += 1
        plotIndex = 1
    

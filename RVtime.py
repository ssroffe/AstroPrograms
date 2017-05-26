import tools as tls
from astropy.io import fits
import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from plot_format import plot_format
import scipy.stats as stats
from datetime import datetime
import time
import os
import sys


def lorentzian(x, depth, width, RVShift):
    return 1.0-depth/(1.0 + ((x-RVShift)/width)**2)

lines = [line.rstrip('\n') for line in open('filelist')]
plot_format()
#print lines
plotNumber = 1
plotIndex = 1

rvJoined = []
numJoined = []
errorJoined = []
timeJoined = []

rvBeta = []
rvGamma = [] 
rvDelta = []

errorBeta = []
errorGamma = []
errorDelta = []

numBeta = []
numGamma = []
numDelta = []

timeBeta = []
timeGamma = []
timeDelta = []

for j in range(len(lines)):

    #plotIndex += 1
    path = lines[j]

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    
    #timeTaken = basename[7:]
    #print timeTaken
    dateTaken = basename[7:17]
    timeTaken = basename[18:23]
    dateTimeTaken = basename[7:23]
    timeFormat = '%Y-%m-%dT%H:%M'

    dateTimeObj = datetime.strptime(dateTimeTaken, timeFormat)
    
    dateTimeNum = dates.date2num(dateTimeObj)
    #print dateTimeNum
    #timeJoined.append(dateTimeObj)
    timeJoined.append(dateTimeNum)

    #print dateTimeObj

    #print "date " + dateTaken
    #print "time " + timeTaken
    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    
    Owl,Normflux,errorNorm = tls.NormNoPlot(path)
    
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
        error = errorNorm[wherr]
        
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

        if lineIndex == 1:
            popt, pcov = sp.curve_fit(lorentzian,vel,flux,maxfev=1500)
            perr = np.sqrt(np.diag(pcov))
            #err = stats.sem(perr)
            err = np.mean(perr)
            rvs = popt[2]

            rvBeta.append(rvs)
            numBeta.append(j)
            errorBeta.append(err)
        if lineIndex == 2:
            popt2, pcov2 = sp.curve_fit(lorentzian,vel,flux,maxfev=1500)
            perr2 = np.sqrt(np.diag(pcov2))
            #err2 = stats.sem(perr2)
            err2 = np.mean(perr2)
            rvs2 = popt2[2]

            rvGamma.append(rvs2)
            numGamma.append(j)
            errorGamma.append(err2)
        if lineIndex == 3:
            popt3, pcov3 = sp.curve_fit(lorentzian,vel,flux,maxfev=1500)
            perr3 = np.sqrt(np.diag(pcov3))
            #err3 = stats.sem(perr3)
            err3 = np.mean(perr3)
            rvs3 = popt3[2]

            rvDelta.append(rvs3)
            numDelta.append(j)
            errorDelta.append(err3)

    errors = np.array([err,err2,err3])
    rvshifts = np.array([rvs,rvs2,rvs3])
    whereMinError = np.argmin(errors)

    
    rvJoined.append(rvshifts[whereMinError])
    numJoined.append(j)
    errorJoined.append(errors[whereMinError])
            
plot_format()
plt.errorbar(timeJoined,rvJoined,yerr=errorJoined,color='black',linestyle='None',marker='o',label="Joined RV shift")

#plt.legend()
plt.ylabel("Velocity [km/s]")
plt.xlabel("Time")
plt.title(wdName+" Joined fits")
plt.gcf().autofmt_xdate()
plt.savefig("RVnumPlot/"+wdName+"_vel_v_time_joined.pdf")
plt.figure()

plot_format()

plt.errorbar(timeJoined,rvBeta,yerr=errorBeta,linestyle='None',marker='o',color='blue',label="H Beta RV shift")
plt.errorbar(timeJoined,rvGamma,yerr=errorGamma,linestyle='None',marker='o',color='green',label="H Gamma RV shift")
plt.errorbar(timeJoined,rvDelta,yerr=errorDelta,linestyle='None',marker='o',color='red',label="H Delta RV shift")

plt.legend()
plt.xlabel("Time")
plt.ylabel("Velocity [km/s]")
plt.title(wdName)
plt.savefig("RVnumPlot/"+wdName+"_vel_v_time_separate.pdf")






##################
#### Old Code ####
##################

"""        
        if rvs == None:
            popt, pcov = sp.curve_fit(lorentzian,vel,flux,maxfev=1500)
            rvs = popt[2]
            perr = np.sqrt(np.diag(pcov))
            rvJoinedArray.append(rvs)
            specNumArray.append(j)
            errorJoinedArray.append(stats.sem(perr))
            
            rvBeta.append(rvs)
            numBeta.append(j)
            errorBeta.append(stats.sem(perr))
            
        else:
            popt, pcov = sp.curve_fit(lambda x,depth,width: lorentzian(x,depth,width,rvs),vel,flux,maxfev=1500)
            perr = np.sqrt(np.diag(pcov))
            popt2, pcov2 = sp.curve_fit(lorentzian,vel,flux,maxfev=1500)
            rvs2 = popt2[2]
            perr2 = np.sqrt(np.diag(pcov2))
"""

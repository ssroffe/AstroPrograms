import tools as tls
from astropy.io import fits
import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from plot_format import plot_format
import os
import sys
import emcee as mc
import corner as cn
from datetime import datetime
from astropy.time import Time

def sine(t,A,P,Phi,Gamma):
    ### A [km/s], P [same as t], Phi [rad], Gamma [km/s]
    return A*np.sin((2*np.pi)*(t/P) + Phi) + Gamma

def lnlikeSine(p,x,y,err):
    A,P,Phi,Gamma = p
    #return -np.sum((y-sine(x,A,P,Phi,Gamma))**2/(2*err))
    return -0.5*np.sum(np.log(err**2)+(y-sine(x,A,P,Phi,Gamma))**2/(err**2))

def lnpriorSine(p):
    A,P,Phi,Gamma = p
    if 0.0 < A < 500.0 and 0.0 < P < 0.1 and 0.0 < Phi < (2*np.pi) and -500.0 < Gamma < 500.0:
        return 0.0
    return -np.inf

def lnprobSine(p,x,y,yerr):
    lp = lnpriorSine(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlikeSine(p,x,y,yerr)

def lorentzian(x, depth, width, RVShift):
    return 1.0-depth/(1.0 + ((x-RVShift)/width)**2)

def lnlike(p,x,y,err):
    depth,width,RVShift = p
    return -np.sum((y-lorentzian(x,depth,width,RVShift))**2/(2*err))

def lnprior(p):
    depth,width,RVShift = p
    if 0.0 < depth < 1.0 and 0.0 < width < 2000 and -500.0 < RVShift < 500.0:
        return 0.0
    return -np.inf

def lnprob(p, x, y, yerr):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p, x, y, yerr)

def lnprobSum(p,x,y,yerr):
    d1,w1,d2,w2,d3,w3,rvs = p
    x1,x2,x3 = x
    y1,y2,y3 = y
    yerr1,yerr2,yerr3 = yerr
    
    p1 = (d1,w1,rvs)
    p2 = (d2,w2,rvs)
    p3 = (d3,w3,rvs)
    s = lnprob(p1,x1,y1,yerr1) + lnprob(p2,x2,y2,yerr2) + lnprob(p3,x3,y3,yerr3)
    return s

lines = [line.rstrip('\n') for line in open('filelist')]
plot_format()
#print lines
plotNumber = 1
plotIndex = 1

numArr = []
rvArr = []
stdArr = []
timeArr = []

for j in range(len(lines)):
#for j in range(0,1):

    path = lines[j]

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    timeTaken = basename[15:]

    dateTaken = basename[7:17]
    timeTaken = basename[18:23]
    dateTimeTaken = basename[7:23]

    t = Time(dateTimeTaken, format='isot',scale='utc')

    timeArr.append(t.mjd)

    ### Get velocites, fluxes and errors
    vels,fluxes,ferrs = tls.GetAllVelocities(path)

    ### Do the fit
    ndim, nwalkers = 7, 100
    sampler = tls.MCMCfit(lnprobSum,args=(vels,fluxes,ferrs),nwalkers=nwalkers,ndim=ndim)

    ### Get the RV Shifts from the sampler with errors
    rvFit,rvStd = tls.GetRV(sampler)

    numArr.append(j)
    rvArr.append(rvFit)
    stdArr.append(rvStd)


timeArr = np.array(timeArr)
numArr = np.array(numArr)
rvArr = np.array(rvArr)
stdArr = np.array(stdArr)

sampler = tls.MCMCfit(lnprobSine,args=(timeArr,rvArr,stdArr),nwalkers=200,ndim=4,burnInSteps=4000,steps=4000)

tls.plotchains(sampler,4,"OrbitalFitsLorentzian/"+wdName,"chains.pdf")

samples = sampler.flatchain.reshape(-1,4).T

AmpArr = samples[0]
PerArr = samples[1]
PhiArr = samples[2]
GamArr = samples[3]

Amp = AmpArr.mean()
Per = PerArr.mean()
Phi = PhiArr.mean()
Gamma = GamArr.mean()

print Amp,Per,Phi,Gamma
#Amp, Per, Phi, Gamma = popt
#newTime = np.linspace(np.min(numArr),np.max(numArr),1000)
newTime = np.linspace(np.min(timeArr),np.max(timeArr),1000)

plot_format()
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="black")
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.title(wdName + " velocity vs time")
plt.savefig("OrbitalFitsLorentzian/"+wdName+"_time.pdf")

plot_format()
plt.errorbar(numArr,rvArr, yerr=stdArr,linestyle='None',marker='o')
#plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="black")
plt.ylabel("Velocity [km/s]")
plt.xlabel("Spectrum number")
plt.title(wdName + " velocity vs spectrum number")
plt.savefig("OrbitalFitsLorentzian/"+wdName+"_num.pdf")

plot_format()
plt.subplot(4,1,1)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="black")
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.title(wdName + " velocity vs time")

plt.subplot(4,1,2)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="black")
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
#plt.title(wdName + " velocity vs time")
plt.xlim(np.min(timeArr)-0.1,np.min(timeArr)+0.1)

plt.subplot(4,1,3)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="black")
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
#plt.title(wdName + " velocity vs time")
plt.xlim(np.median(timeArr)-0.1,np.median(timeArr)+0.1)

plt.subplot(4,1,4)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="black")
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
#plt.title(wdName + " velocity vs time")
plt.xlim(np.max(timeArr)-0.1,np.max(timeArr)+0.1)

plt.savefig("OrbitalFitsLorentzian/"+wdName+"_time_zoomed.pdf")


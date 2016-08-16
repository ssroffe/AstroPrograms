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

def voigt(x, Ldepth, Lwidth, Gdepth, Gwidth, RVShift):
    return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))

def lorentzian(x, depth, width, RVShift):
    return 1.0-depth/(1.0 + ((x-RVShift)/width)**2)

def lnlike(p,x,y,err):
    Ldepth,Lwidth,Gdepth,Gwidth,RVShift = p
    return -np.sum((y-voigt(x,Ldepth,Lwidth,Gdepth,Gwidth,RVShift))**2/(2*err))

def lnprior(p):
    Ldepth,Lwidth,Gdepth,Gwidth,RVShift = p
    if 0.0 < Ldepth < 1.0 and 0.0 < Lwidth < 2000.0 and 0.0 < Gdepth < 1.0 and 0.0 < Gwidth < 2000.0 and -500.0 < RVShift < 500.0:
        return 0.0
    return -np.inf

def lnprob(p, x, y, yerr):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p, x, y, yerr)

def lnprobSum(p, x, y, yerr):
    ld1,wd1,gd1,gw1,ld2,wd2,gd2,gw2,ld3,wd3,gd3,gw3,rvs = p
    x1,x2,x3 = x
    y1,y2,y3 = y
    yerr1,yerr2,yerr3 = yerr
    
    p1 = (ld1,wd1,gd1,gw1,rvs)
    p2 = (ld2,wd2,gd2,gw2,rvs)
    p3 = (ld3,wd3,gd3,gw3,rvs)
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
    #vels,fluxes,ferrs = tls.GetAllVelocities(path)
    vels,fluxes,ferrs = tls.GetVelArr(path,1)

    ### Do the fit
    popt, pcov = sp.curve_fit(lorentzian,vels,fluxes,sigma=ferrs,maxfev=8000)
    perr = np.sqrt(np.diag(pcov))

    rvFit = popt[-1]
    rvStd = perr[-1]

    print rvFit
    print rvStd
    ### Get the RV Shifts from the sampler with errors
    #rvFit,rvStd = tls.GetRV(sampler)

    numArr.append(j)
    rvArr.append(rvFit)
    stdArr.append(rvStd)



timeArr = np.array(timeArr)
numArr = np.array(numArr)
rvArr = np.array(rvArr)
stdArr = np.array(stdArr)

walkers, dim = 200, 4

sampler = tls.MCMCfit(lnprobSine,args=(timeArr,rvArr,stdArr),nwalkers=walkers,ndim=dim,burnInSteps=8000,steps=8000)

tls.plotchains(sampler,4,"OrbitalFitsScipy/"+wdName,"chains.pdf")

samplesChain = sampler.chain[:,:,:].reshape((-1,4))
samples = sampler.flatchain.reshape(-1,4).T

AmpArr = samples[0]
PerArr = samples[1]
PhiArr = samples[2]
GamArr = samples[3]

Amp = AmpArr.mean()
Per = PerArr.mean()
Phi = PhiArr.mean()
Gamma = GamArr.mean()

newTime = np.linspace(np.min(timeArr),np.max(timeArr),10000)

AArr = []
PArr = []
PhArr = []
GArr = []

#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=1000)]:
#    AArr.append(A)
#    PArr.append(P)
#    PhArr.append(Ph)
#    GArr.append(Gam)
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.01)
    #print A,P,Ph,Gam

plot_format()
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
#plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.ylim(-400,400)
#plt.title(wdName + " velocity vs time\n{v[t]="+str(params[0])+"*sin[2*pi*(t/"+str(params[1])+") + "+str(params[2])+"] + "+str(params[3])+"}")
plt.savefig("OrbitalFitsScipy/"+wdName+"_time.pdf")

plot_format()
plt.subplot(4,1,1)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=1000)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.01)
    AArr.append(A)
    PArr.append(P)
    PhArr.append(Ph)
    GArr.append(Gam)
#    print A,P,Ph,Gam

params = [AArr[-1],PArr[-1],PhArr[-1],GArr[-1]]

plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
#plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.ylim(-400,400)
plt.title(wdName + " velocity vs time"+'\n'+"$v(t) = {0:.3f}*sin(2\pi(t/{1:.3f}) + {2:.3f}) + {3:.3f}$".format(params[0],params[1],params[2],params[3]))

plt.subplot(4,1,2)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.ylim(-400,400)
#plt.title(wdName + " velocity vs time")
plt.xlim(np.min(timeArr)-0.1,np.min(timeArr)+0.1)

plt.subplot(4,1,3)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.ylim(-400,400)
#plt.title(wdName + " velocity vs time")
plt.xlim(np.median(timeArr)-0.1,np.median(timeArr)+0.1)

plt.subplot(4,1,4)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.ylim(-400,400)
#plt.title(wdName + " velocity vs time")
plt.xlim(np.max(timeArr)-0.1,np.max(timeArr)+0.1)

plt.savefig("OrbitalFitsScipy/"+wdName+"_time_zoomed.pdf")

#print ""
#print Amp,Per,Phi,Gamma

import corner
fig = corner.corner(samplesChain, labels=["A","P","Phi","Gam"])
fig.savefig("OrbitalFitsScipy/"+wdName+"_Triangle.pdf")


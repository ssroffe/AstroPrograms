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
    return -np.sum((y-sine(x,A,P,Phi,Gamma))**2/(2*err))
    #return -0.5*np.sum(np.log(err**2)+(y-sine(x,A,P,Phi,Gamma))**2/(err**2))

def lnpriorSine(p):
    A,P,Phi,Gamma = p
    #minimum at 0.02 days for P
    if 5.0 < A < 500.0 and 0.02 < P < 0.04 and 0.0 < Phi < (2*np.pi) and -500.0 < Gamma < 500.0:
        return 0.0
    return -np.inf

def lnprobSine(p,x,y,yerr):
    lp = lnpriorSine(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlikeSine(p,x,y,yerr)

def voigt(x, Ldepth, Lwidth, Gdepth, Gwidth, RVShift):
    return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))

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
    ld1,lw1,gd1,gw1,ld2,lw2,gd2,gw2,ld3,lw3,gd3,gw3,rvs = p
    x1,x2,x3 = x
    y1,y2,y3 = y
    yerr1,yerr2,yerr3 = yerr
    
    p1 = (ld1,lw1,gd1,gw1,rvs)
    p2 = (ld2,lw2,gd2,gw2,rvs)
    p3 = (ld3,lw3,gd3,gw3,rvs)
    s = lnprob(p1,x1,y1,yerr1) + lnprob(p2,x2,y2,yerr2) + lnprob(p3,x3,y3,yerr3)
    return s

tls.mkdir("ParamVsNumPlot")

lines = [line.rstrip('\n') for line in open('filelist')]
plot_format()
#print lines
plotNumber = 1
plotIndex = 1

numArr = []
rvArr = []
stdArr = []
timeArr = []

ld1Arr = []
lw1Arr = []
gd1Arr = []
gw1Arr = []

ld2Arr = []
lw2Arr = []
gd2Arr = []
gw2Arr = []

ld3Arr = []
lw3Arr = []
gd3Arr = []
gw3Arr = []


for j in range(len(lines)):
#for j in range(0,1):

    path = lines[j]

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    #timeTaken = basename[15:]

    dateTaken = basename[7:17]
    #timeTaken = basename[18:23]
    #dateTimeTaken = basename[7:23]
    dateTimeTaken = tls.GetDateTime(path)
    if "T" not in dateTimeTaken:
        dateTimeTaken = dateTaken + "T" + dateTimeTaken

    t = Time(dateTimeTaken, format='isot',scale='utc')
    
    timeArr.append(t.mjd)

    ### Get velocites, fluxes and errors
    vels,fluxes,ferrs = tls.GetAllVelocities(path)

    ### Do the fit
    #ndim, nwalkers = 7,200
    ndim, nwalkers = 13, 200
    
    sampler = tls.MCMCfit(lnprobSum,args=(vels,fluxes,ferrs),nwalkers=nwalkers,ndim=ndim,burnInSteps=4000,steps=4000)

    samples = sampler.flatchain.reshape(-1,ndim).T

    ld1Arr.append(samples[0].mean())
    lw1Arr.append(samples[1].mean())
    gd1Arr.append(samples[2].mean())
    gw1Arr.append(samples[3].mean())

    ld2Arr.append(samples[4].mean())
    lw2Arr.append(samples[5].mean())
    gd2Arr.append(samples[6].mean())
    gw2Arr.append(samples[7].mean())

    ld3Arr.append(samples[8].mean())
    lw3Arr.append(samples[9].mean())
    gd3Arr.append(samples[10].mean())
    gw3Arr.append(samples[11].mean())
    
    ### Get the RV Shifts from the sampler with errors
    rvFit,rvStd = tls.GetRV(sampler)
        

    numArr.append(j)
    rvArr.append(rvFit)
    stdArr.append(rvStd)
    

timeArr = np.array(timeArr)
numArr = np.array(numArr)
rvArr = np.array(rvArr)
stdArr = np.array(stdArr)

plot_format()
plt.plot(numArr,ld1Arr,color="blue",marker="o",label="L-depth [HBeta]",linestyle="None")
plt.plot(numArr,ld2Arr,color="red",marker="o",label="L-depth [HGamma]",linestyle="None")
plt.plot(numArr,ld3Arr,color="green",marker="o",label="L-depth [HDelta]",linestyle="None")
plt.plot(numArr,gd1Arr,color="blue",marker="s",label="G-depth [HBeta]",linestyle="None")
plt.plot(numArr,gd2Arr,color="red",marker="s",label="G-depth [HGamma]",linestyle="None")
plt.plot(numArr,gd3Arr,color="green",marker="s",label="G-depth [HDelta]",linestyle="None")

plt.xlabel("Spectrum Number")
plt.ylabel("Depth")
plt.title(wdName+" Depth Vs Spectrum Number")
plt.legend()
plt.savefig("ParamVsNumPlot/"+wdName+"_DepthVsNum.pdf")

plot_format()

plt.plot(numArr,lw1Arr,color="blue",marker="o",label="L-width [HBeta]",linestyle="None")
plt.plot(numArr,lw2Arr,color="red",marker="o",label="L-width [HGamma]",linestyle="None")
plt.plot(numArr,lw3Arr,color="green",marker="o",label="L-width [HDelta]",linestyle="None")
plt.plot(numArr,gw1Arr,color="blue",marker="s",label="G-width [HBeta]",linestyle="None")
plt.plot(numArr,gw2Arr,color="red",marker="s",label="G-width [HGamma]",linestyle="None")
plt.plot(numArr,gw3Arr,color="green",marker="s",label="G-width [HDelta]",linestyle="None")

plt.xlabel("Spectrum Number")
plt.ylabel("Width")
plt.title(wdName+" Width Vs Spectrum Number")
plt.legend()

plt.savefig("ParamVsNumPlot/"+wdName+"_WidthVsNum.pdf")


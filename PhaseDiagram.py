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
    ld1,wd1,gd1,gw1,ld2,wd2,gd2,gw2,ld3,wd3,gd3,gw3,rvs = p
    x1,x2,x3 = x
    y1,y2,y3 = y
    yerr1,yerr2,yerr3 = yerr
    
    p1 = (ld1,wd1,gd1,gw1,rvs)
    p2 = (ld2,wd2,gd2,gw2,rvs)
    p3 = (ld3,wd3,gd3,gw3,rvs)
    s = lnprob(p1,x1,y1,yerr1) + lnprob(p2,x2,y2,yerr2) + lnprob(p3,x3,y3,yerr3)
    return s


tls.mkdir("PhaseDiagram")

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


    ### Get the RV Shifts from the sampler with errors
    rvFit,rvStd = tls.GetRV(sampler)
        

    numArr.append(j)
    rvArr.append(rvFit)
    stdArr.append(rvStd)
    

timeArr = np.array(timeArr)
numArr = np.array(numArr)
rvArr = np.array(rvArr)
stdArr = np.array(stdArr)

Amin, Amax = 5.0, 500.0
Pmin, Pmax = 0.02,0.04
Phimin,Phimax = 0.0,(2*np.pi)
GamMin, GamMax = -500.0, 500.0

middles = np.array([(Amin+Amax)/2,(Pmin+Pmax)/2,(Phimin+Phimax)/2,(GamMin+GamMax)/2])

walkers,dim = 200,4

pos = [middles + 1e-4*np.random.randn(dim) for i in range(walkers)]


sampler = tls.MCMCfit(lnprobSine,args=(timeArr,rvArr,stdArr),nwalkers=walkers,ndim=4,burnInSteps=250000,steps=250000,p=pos)

samplesChain = sampler.chain[:,:,:].reshape((-1,4))
samples = sampler.flatchain.reshape((-1,4)).T

AArr = []
PArr = []
PhArr = []
GArr = []
for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=1000)]:
    AArr.append(A)
    PArr.append(P)
    PhArr.append(Ph)
    GArr.append(Gam)
    
params = [AArr[-1],PArr[-1],PhArr[-1],GArr[-1]]

AFit,PFit,PhFit,GamFit = params

phiDiagArr = []

#A*np.sin((2*np.pi)*(t/P) + Phi) + Gamma
for pt in timeArr:
    PhiOff = ((2*np.pi)*(pt/PFit) + PhFit) 
    phiDiagArr.append(PhiOff)

### make sure the points are 0 < point < 2pi
for i in range(len(phiDiagArr)):
    while (phiDiagArr[i] < 0) or (phiDiagArr[i] > (2*np.pi)):
        if (phiDiagArr[i] < 0):
            phiDiagArr[i] = phiDiagArr[i] + (2*np.pi)
        else:
            phiDiagArr[i] = phiDiagArr[i] - (2*np.pi)
    
phiDiagArr = np.array(phiDiagArr)

angles = np.linspace(0,2*np.pi,5000)

yvalues = sine(phiDiagArr,AFit,2*np.pi,0,GamFit)

residuals = yvalues - rvArr


#### PLOTTING STUFF

plt.figure(1).add_axes((.1,.3,.8,.6))
plt.plot(angles,sine(angles,AFit,2*np.pi,0,GamFit),'k--')
plt.title(wdName+" Phase")
plt.ylabel("RV [km/s]")
plt.xlim(0,2*np.pi)
plt.errorbar(phiDiagArr,rvArr,yerr=stdArr,linestyle='None',marker='o')

plt.figure(1).add_axes((.1,.1,.8,.2))
plt.plot(phiDiagArr,residuals,linestyle="None",marker="o")
plt.xlabel("angle [rad]")
plt.xlim(0,2*np.pi)
plt.ylim(-200,200)
plt.axhline(0,linestyle="--",color="black",alpha=0.5)
plt.savefig("PhaseDiagram/"+wdName+"_phase.pdf")

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
from paperplots import Signal2Noise

def NoOrbit(t,Gamma):
    return Gamma

def lnlikeNoOrbit(p,x,y,err):
    Gamma = p
    return -np.sum((y-NoOrbit(x,Gamma))**2/(2*err))

def lnpriorNoOrbit(p):
    Gamma = p
    if -500.0 < Gamma < 500.0:
        return 0.0
    return -np.inf

def lnprobNoOrbit(p,x,y,yerr):
    lp = lnpriorNoOrbit(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlikeNoOrbit(p,x,y,yerr)

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
    if 5.0 < A < 500.0 and 0.02 < P < 0.1 and 0.0 < Phi < (2*np.pi) and -500.0 < Gamma < 500.0:
        return 0.0
    return -np.inf

def lnprobSine(p,x,y,yerr):
    lp = lnpriorSine(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlikeSine(p,x,y,yerr)

if (len(sys.argv) == 2) and sys.argv[1] is "lorentzian":
    print "Using Lorentzian model"
    def lorentzianModel(x, Ldepth, Lwidth):
        return 1.0-Ldepth/(1.0 + ((x)/Lwidth)**2)

    def lnlikeModel(p,x,y,err):
        Ldepth,Lwidth = p
        return -np.sum((y-lorentzianModel(x,Ldepth,Lwidth))**2/(2*err))

    def lnpriorModel(p):
        Ldepth,Lwidth = p
        if 0.0 < Ldepth < 1.0 and 0.0 < Lwidth < 2000.0:
            return 0.0
        return -np.inf

    def lnprobModel(p, x, y,yerr):
        lp = lnpriorModel(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlikeModel(p, x, y, yerr)


    def lorentzian(x, Ldepth, Lwidth, RVShift):
        return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2)

    def lnlike(p,x,y,err):
        Ldepth,Lwidth,RVShift = p
        return -np.sum((y-lorentzian(x,Ldepth,Lwidth,RVShift))**2/(2*err))

    def lnprior(p):
        Ldepth,Lwidth,RVShift = p
        if -500.0 < RVShift < 500.0:
            return 0.0
        return -np.inf

    def lnprob(p, x, y, yerr):
        lp = lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(p, x, y, yerr)

    def lnprobSum(p, x, y, yerr):
        rvs = p
        
        x1,x2,x3 = x
        y1,y2,y3 = y
        yerr1,yerr2,yerr3 = yerr
        
        p1 = (ld1,lw1,rvs)
        p2 = (ld2,lw2,rvs)
        p3 = (ld3,lw3,rvs)
        s = lnprob(p1,x1,y1,yerr1) + lnprob(p2,x2,y2,yerr2) + lnprob(p3,x3,y3,yerr3)
        return s

else:
    print "Using Voigt Model"
    def voigtModel(x, Ldepth, Lwidth, Gdepth, Gwidth):
        return 1.0-Ldepth/(1.0 + ((x)/Lwidth)**2) - Gdepth*np.exp(-(x)**2/(2*Gwidth**2))

    def lnlikeModel(p,x,y,err):
        Ldepth,Lwidth,Gdepth,Gwidth = p
        return -np.sum((y-voigtModel(x,Ldepth,Lwidth,Gdepth,Gwidth))**2/(2*err))

    def lnpriorModel(p):
        Ldepth,Lwidth,Gdepth,Gwidth = p
        if 0.0 < Ldepth < 1.0 and 0.0 < Lwidth < 3000.0 and 0.0 < Gdepth < 1.0 and 0.0 < Gwidth < 600.0:
            return 0.0
        return -np.inf
    
    def lnprobModel(p, x, y,yerr):
        lp = lnpriorModel(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlikeModel(p, x, y, yerr)

    
    def voigt(x, Ldepth, Lwidth, Gdepth, Gwidth, RVShift):
        return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))
    
    def lnlike(p,x,y,err):
        Ldepth,Lwidth,Gdepth,Gwidth,RVShift = p
        return -np.sum((y-voigt(x,Ldepth,Lwidth,Gdepth,Gwidth,RVShift))**2/(2*err))
    
    def lnprior(p):
        Ldepth,Lwidth,Gdepth,Gwidth,RVShift = p
        if -500.0 < RVShift < 500.0:
            return 0.0
        return -np.inf

    def lnprob(p, x, y, yerr):
        lp = lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(p, x, y, yerr)

    def lnprobSum(p, x, y, yerr):
        rvs = p
        
        x1,x2,x3 = x
        y1,y2,y3 = y
        yerr1,yerr2,yerr3 = yerr
        """print "\n"
        print ld1, lw1, gd1, gw1
        print ld2, lw2, gd2, gw2
        print ld3, lw3, gd3, gw3
        print "\n"
        """
        p1 = (ld1,lw1,gd1,gw1,rvs)
        p2 = (ld2,lw2,gd2,gw2,rvs)
        p3 = (ld3,lw3,gd3,gw3,rvs)
        s = lnprob(p1,x1,y1,yerr1) + lnprob(p2,x2,y2,yerr2) + lnprob(p3,x3,y3,yerr3)
        
        return s

#tls.mkdir("BICFits")
#tls.mkdir("BICFits/VelFits")

#tls.mkdir("AICFits")

lines = [line.rstrip('\n') for line in open('filelist')]

plot_format()

modelFile = "../../KoesterModels/da"+open('modelVal').read().splitlines()[0]+".dk"

modelWl,modelFlux = tls.ModelNormNoPlot(modelFile)
tmpWl,tmpFlux,_ = tls.NormNoPlot(lines[0])
plot_format()
plt.plot(tmpWl,tmpFlux,alpha=0.4)
plt.plot(modelWl,modelFlux)
plt.xlim(np.min(tmpWl),np.max(tmpWl))
plt.savefig("AICFits/SpectrumCheck.pdf")
plot_format()

modelVels,modelFluxes = tls.ModelGetAllVelocities(modelFile)

modelErrs = []
for i in range(len(modelFluxes)):
    modelErrs.append(0.01*modelFluxes[i])

plotNumber = 1
plotIndex = 1

if (len(sys.argv) == 2) and sys.argv[1] is "lorentzian":
    mdim,mwalkers = 2,200
else:
    mdim,mwalkers = 4,200
"""
for i in range(len(modelVels)):
    
    modelSampler = tls.MCMCfit(lnprobModel,args=(np.array(modelVels[i]),np.array(modelFluxes[i]),np.array(modelErrs[i])),nwalkers=mwalkers,ndim=mdim,burnInSteps=16000,steps=16000)
    modelSamples = modelSampler.flatchain.reshape((-1,mdim)).T
    if i == 0:
        global ld1
        ld1 = modelSamples[0].mean()
        global lw1
        lw1 = modelSamples[1].mean()
        if (len(sys.argv) == 1) or (sys.argv[1] is "voigt"):
            global gd1
            gd1 = modelSamples[2].mean()
            global gw1
            gw1 = modelSamples[3].mean()
    elif i == 1:
        global ld2
        ld2 = modelSamples[0].mean()
        global lw2
        lw2 = modelSamples[1].mean()
        if (len(sys.argv) == 1) or (sys.argv[1] is "voigt"):
            global gd2
            gd2 = modelSamples[2].mean()
            global gw2
            gw2 = modelSamples[3].mean()
    elif i == 2:
        global ld3
        ld3 = modelSamples[0].mean()
        global lw3
        lw3 = modelSamples[1].mean()
        if (len(sys.argv) == 1) or (sys.argv[1] is "voigt"):
            global gd3
            gd3 = modelSamples[2].mean()
            global gw3
            gw3 = modelSamples[3].mean() 


numArr = []
rvArr = []
stdArr = []
timeArr = []
"""
plot_format()
for j in range(len(lines)):
#for j in range(0,1):
    
    path = lines[j]

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    #timeTaken = basename[15:]

rvdata = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')

timeArr = rvdata[:,0]
rvArr = rvdata[:,1]
stdArr = rvdata[:,2]

SNArr = Signal2Noise()

SNCut = 2.0
wherrSN = np.where(SNArr >= SNCut)

timeArr = timeArr[wherrSN]
rvArr = rvArr[wherrSN]
stdArr = stdArr[wherrSN]

###
Amin, Amax = 5.0, 500.0
Pmin, Pmax = 0.02,0.1
Phimin,Phimax = 0.0,(2*np.pi)
GamMin, GamMax = -500.0, 500.0

#plot_format()
#plt.errorbar(numArr,rvArr,yerr=stdArr,ls='None')
#plt.ylabel("RV [km/s]")
#plt.xlabel("Spectrum number")
#plt.savefig("BICFits/"+wdName+"_numRV.pdf")
#plot_format()

middles = np.array([(Amin+Amax)/2,(Pmin+Pmax)/2,(Phimin+Phimax)/2,(GamMin+GamMax)/2])

#noOrbWalkers,noOrbDim = 200,1
#noOrbPos = [middles[-1] + 1e-4*np.random.randn(noOrbDim) for i in range(noOrbWalkers)]

#noOrbSampler = tls.MCMCfit(lnprobNoOrbit,args=(timeArr,rvArr,stdArr),nwalkers=noOrbWalkers,ndim=noOrbDim,burnInSteps=250000,steps=250000,p=noOrbPos)

#noOrbSamplesChain = noOrbSampler.chain[:,:,:].reshape((-1,1))
#noOrbSamples = noOrbSampler.flatchain.reshape((-1,1)).T

#mArr = []
#for m in noOrbSamplesChain[np.random.randint(len(noOrbSamplesChain),size=1000)]:
#    mArr.append(m)
#minit = mArr[-1]

#nllNoOrb = lambda *args: -lnprobNoOrbit(*args)
#mparam = sp.minimize(nllNoOrb,[minit],args=(timeArr,rvArr,stdArr))["x"]
#print mparam

wgtAvg = (np.sum(rvArr * stdArr**(-2))) / np.sum(stdArr[i]**(-2))
wgtStd = 1 / (np.sum(stdArr**(-2)))
mparam = wgtAvg

walkers,dim = 200,4

pos = [middles + 1e-4*np.random.randn(dim) for i in range(walkers)]

sampler = tls.MCMCfit(lnprobSine,args=(timeArr,rvArr,stdArr),nwalkers=walkers,ndim=4,burnInSteps=250000,steps=250000,p=pos)

tls.plotchains(sampler,4,"OrbitalFitsVoigt/"+wdName,"chains.pdf")

samplesChain = sampler.chain[:,:,:].reshape((-1,4))
samples = sampler.flatchain.reshape((-1,4)).T

#print np.shape(samplesChain)
#print samples

AmpArr = samples[0]
PerArr = samples[1]
PhiArr = samples[2]
GamArr = samples[3]

Amp = AmpArr.mean()
Per = PerArr.mean()
Phi = PhiArr.mean()
Gamma = GamArr.mean()

PhiStd = PhiArr.std()

newTime = np.linspace(np.min(timeArr),np.max(timeArr),5000)

NoOrbArr = []
for i in range(len(newTime)):
    NoOrbArr.append(mparam)
NoOrbArr = np.array(NoOrbArr)

AArr = []
PArr = []
PhArr = []
GArr = []

plot_format()
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
#plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
plt.plot(newTime,NoOrbArr,color='r')
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.savefig("AICFits/"+wdName+"_time.pdf")

plot_format()
plt.subplot(4,1,1)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=5000)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.01)
    AArr.append(A)
    PArr.append(P)
    PhArr.append(Ph)
    GArr.append(Gam)
    #print A,P,Ph,Gam

nll = lambda *args: -lnprobSine(*args)
results = sp.minimize(nll, [AArr[-1],PArr[-1],PhArr[-1],GArr[-1]],args=(timeArr,rvArr,stdArr))

Astd = np.array(AArr).std()
Pstd = np.array(PArr).std()
Phstd = np.array(PhArr).std()
Gstd = np.array(GArr).std()

#print results
params = []
Afit,Pfit,Phfit,Gfit = results["x"]
#params = [AArr[-1],PArr[-1],PhArr[-1],GArr[-1]]
params = [(Afit,Astd),(Pfit,Pstd),(Phfit,Phstd),(Gfit,Gstd)]

np.savetxt("AICFits/"+wdName+"_sineParams.csv",params,delimiter=',')

##### BIC CALCULATIONS ########
noOrbParams = (mparam)
noOrbk = 1
#noOrbBIC = -2*lnlikeNoOrbit(mparam,timeArr,rvArr,stdArr)+noOrbk*np.log(len(timeArr))
noOrbBIC = -2*lnlikeNoOrbit(mparam,timeArr,rvArr,stdArr)+2*noOrbk + ( (2*noOrbk*(noOrbk+1)) / (len(timeArr) - noOrbk - 1))

sineParams = (Afit,Pfit,Phfit,Gfit)
sinek = 4
#sineBIC = -2*lnlikeSine(sineParams,timeArr,rvArr,stdArr)+sinek*np.log(len(timeArr))
sineBIC = -2*lnlikeNoOrbit(mparam,timeArr,rvArr,stdArr)+2*sinek + ( (2*sinek*(sinek+1)) / (len(timeArr) - sinek - 1))

deltaBIC = noOrbBIC - sineBIC

bicFile = open("AICFits/"+wdName+"_BICCalc.txt",'w')
bicFile.write("Orbit eqn: v(t) = {0:.3f}*sin(2*pi*(t/{1:.3f}) + {2:.3f}) + {3:.3f}\n".format(Afit,Pfit,Phfit,Gfit))
bicFile.write("No Orbit eqn: v(t) = {0:.3f}\n".format(float(mparam)))
bicFile.write("No Orbit BIC = {0:.3f}\n".format(float(noOrbBIC)))
bicFile.write("Sine AIC = {0:.3f}\n".format(float(sineBIC)))
bicFile.write("Delta AIC = noOrbBIC - sineBIC = {0:.3f}".format(float(deltaBIC)))
bicFile.close()
deltaBICArr = np.array([noOrbBIC, sineBIC, deltaBIC])
np.savetxt("AICFits/"+wdName+"_deltaBIC.csv",deltaBICArr,delimiter=',')
##############################


plt.plot(newTime,sine(newTime,Afit,Pfit,Phfit,Gfit),color="red",alpha=0.75)
#plt.plot(newTime,sine(newTime,Amp,Per,Phi,Gamma),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
plt.title(wdName + " velocity vs time"+'\n'+"$v(t) = {0:.3f}*sin(2\pi(t/{1:.3f}) + {2:.3f}) + {3:.3f}$".format(Afit,Pfit,Phfit,Gfit))

plt.subplot(4,1,2)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,Afit,Pfit,Phfit,Gfit),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
#plt.title(wdName + " velocity vs time")
plt.xlim(np.min(timeArr)-0.1,np.min(timeArr)+0.1)

plt.subplot(4,1,3)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,Afit,Pfit,Phfit,Gfit),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
#plt.title(wdName + " velocity vs time")
plt.xlim(np.median(timeArr)-0.1,np.median(timeArr)+0.1)

plt.subplot(4,1,4)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,Afit,Pfit,Phfit,Gfit),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
#plt.title(wdName + " velocity vs time")
plt.xlim(np.max(timeArr)-0.1,np.max(timeArr)+0.1)

plt.savefig("AICFits/"+wdName+"_time_zoomed.pdf")

#print ""
#print Amp,Per,Phi,Gamma

import corner
fig = corner.corner(samplesChain, labels=["A","P","Phi","Gam"])
fig.savefig("AICFits/"+wdName+"_Triangle.pdf")


#!/usr/bin/python

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

import corner

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
        if 0.0 < Ldepth < 1.0 and 1000.0 < Lwidth < 3000.0 and 0.0 < Gdepth < 1.0 and 10.0 < Gwidth < 300.0:
            return 0.0
        return -np.inf
    
    def lnprobModel(p, x, y,yerr):
        lp = lnpriorModel(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlikeModel(p, x, y, yerr)

    
    def voigt(x, Ldepth, Lwidth, Gdepth, Gwidth, RVShift):
        return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))
    def voigt1(x, RVShift):
        Ldepth,Lwidth,Gdepth,Gwidth = ld1,lw1,gd1,gw1
        return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))
    def voigt2(x, RVShift):
        Ldepth,Lwidth,Gdepth,Gwidth = ld2,lw2,gd2,gw2
        return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))
    def voigt3(x, RVShift):
        Ldepth,Lwidth,Gdepth,Gwidth = ld3,lw3,gd3,gw3
        return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))
    
    #def lnlike(p,x,y,err):
    #    Ldepth,Lwidth,Gdepth,Gwidth,RVShift = p
    #    return -np.sum((y-voigt(x,Ldepth,Lwidth,Gdepth,Gwidth,RVShift))**2/(2*err))

    def lnlike1(p,x,y,err):
        RVShift = p
        return -np.sum((y-voigt1(x,RVShift))**2/(2*err))
    def lnlike2(p,x,y,err):
        RVShift = p
        return -np.sum((y-voigt2(x,RVShift))**2/(2*err))
    def lnlike3(p,x,y,err):
        RVShift = p
        return -np.sum((y-voigt3(x,RVShift))**2/(2*err))
    
    
    def lnprior(p):
        RVShift = p
        #if -500.0 < RVShift < 500.0:
        if tmpsdssRV == None and tmpsdssSTD == None:
            if -250.0 < RVShift < 250.0:
                return 0.0
            return -np.inf
        else:
            if (tmpsdssRV - tmpsdssSTD) < RVShift < (tmpsdssRV - tmpsdssSTD):
                return 0.0
            return -np.inf

    def lnprob(p, x, y, yerr):
        lp = lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(p, x, y, yerr)
    def lnprob1(p, x, y, yerr):
        lp = lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike1(p, x, y, yerr)
    def lnprob2(p, x, y, yerr):
        lp = lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike2(p, x, y, yerr)
    def lnprob3(p, x, y, yerr):
        lp = lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike3(p, x, y, yerr)

    
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
        #s = lnprob(p1,x1,y1,yerr1) + lnprob(p2,x2,y2,yerr2) + lnprob(p3,x3,y3,yerr3)
        s = lnprob1(p,x1,y1,yerr1) + lnprob2(p,x2,y2,yerr2) + lnprob3(p,x3,y3,yerr3)
        
        return s


    def lnprior2(p):
        RVShift = p
        #if -500.0 < RVShift < 500.0:
        if (tmpsdssRV-tmpsdssSTD) < RVShift < (tmpsdssRV+tmpsdssSTD):
            return 0.0
        return -np.inf

global tmpsdssRV, tmpsdssSTD
tmpsdssRV,tmpsdssSTD = None,None
lines = [line.rstrip('\n') for line in open('filelist')]

plot_format()

modelFile = "../../KoesterModels/da"+open('modelVal').read().splitlines()[0]+".dk"

sdssFile = "tmpSpec"

sdssData = np.genfromtxt(sdssFile,skip_header=1,delimiter=',')

sdssWl = sdssData[:,0]
sdssFlux = np.array(sdssData[:,1])

sdssFlux = sdssFlux / max(sdssFlux)

modelData = np.genfromtxt(modelFile,skip_header=34)
modelFlux = np.array(modelData[:,1])
modelWl = modelData[:,0]

modelFlux = modelFlux / max(modelFlux)

tls.mkdir("sdssFits")
plot_format()
plt.plot(sdssWl,sdssFlux,label="sdss")
plt.plot(modelWl,modelFlux,label="Model")
plt.legend()
plt.savefig("sdssFits/sdss_v_Model_raw.pdf")
#plt.show()
plot_format()

modelWl,modelFlux = tls.ModelNormNoPlot(modelFile)
tmpWl,tmpFlux,_ = tls.NormNoPlot(lines[0])

sdssWl, sdssFlux = tls.CSVNormNoPlot(sdssFile)

plt.plot(sdssWl,sdssFlux,alpha=0.5)
plt.plot(modelWl,modelFlux)
plt.xlim(min(tmpWl),max(tmpWl))
plt.axhline(1,color='k',ls='--')
plt.ylim(0,2)
plt.savefig("sdssFits/sdss_v_Model_normalized.pdf")
#plt.show()

sdssVels, sdssFluxes = tls.CSVGetAllVelocities(sdssFile)

sdssErrs = []
for i in range(len(sdssFluxes)):
    sdssErrs.append(0.1*sdssFluxes[i])


    
plot_format()
plt.plot(tmpWl,tmpFlux,alpha=0.4)
plt.plot(modelWl,modelFlux)
plt.xlim(np.min(tmpWl),np.max(tmpWl))
plt.savefig("ModelFits/SpectrumCheck.pdf")
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

LDMin,LDMax = 0.0,1.0
LWMin,LWMax = 1000.0,3000.0
GDMin,GDMax = 0.0,1.0
GWMin,GWMax = 10.0, 300.0
    
modelMiddles = np.array([(LDMin + LDMax)/2, (LWMin + LWMax)/2, (GDMin + GDMax)/2,(GWMin + GWMax)/2])
    
modelPos = [modelMiddles + 1e-4*np.random.randn(mdim) for i in range(mwalkers)]
    
for i in range(len(modelVels)):
    
    modelSampler = tls.MCMCfit(lnprobModel,args=(np.array(modelVels[i]),np.array(modelFluxes[i]),np.array(modelErrs[i])),nwalkers=mwalkers,ndim=mdim,burnInSteps=16000,steps=16000,p=modelPos)
    modelSamplesChain = modelSampler.chain[:,:,:].reshape((-1,mdim))
    modelSamples = modelSampler.flatchain.reshape((-1,mdim)).T
    plot_format()
    plt.plot(modelVels[i],modelFluxes[i],alpha=0.5)
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
        plt.plot(modelVels[i],voigtModel(modelVels[i],ld1,lw1,gd1,gw1))

        #modelSampler.reset()
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
        plt.plot(modelVels[i],voigtModel(modelVels[i],ld2,lw2,gd2,gw2))
        #modelSampler.reset()
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
        plt.plot(modelVels[i],voigtModel(modelVels[i],ld3,lw3,gd3,gw3))

    plt.savefig("sdssFits/ModelLine_"+str(i)+".pdf")
    cor = corner.corner(modelSamplesChain,labels=["LD","LW","GD","GW"])
    cor.savefig("sdssFits/corner_"+str(i)+".pdf")
    modelSampler.reset()


        
tls.mkdir("ModelFits")
tls.mkdir("ModelFits/VelFits")

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
    ndim, nwalkers = 1, 200

    tmpsdssSampler = tls.MCMCfit(lnprobSum,args=(np.array(sdssVels),np.array(sdssFluxes),np.array(sdssErrs)),nwalkers=nwalkers,ndim=ndim,burnInSteps=16000,steps=16000)
    #sdssSamplesChain = sdssSampler.chain[:,:,:].reshape((-1,ndim))
    
    off = 0.5
    tmpsdssRV, tmpsdssSTD = tls.GetRV(tmpsdssSampler)

    tmpMid = np.array([tmpsdssRV])
    sdssPos = [tmpMid + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    sdssSampler = tls.MCMCfit(lnprobSum,args=(np.array(sdssVels),np.array(sdssFluxes),np.array(sdssErrs)),nwalkers=nwalkers,ndim=ndim,burnInSteps=16000,steps=16000,p=sdssPos)
    
    sdssSamplesChain = sdssSampler.chain[:,:,:].reshape((-1,ndim))

    sdssRV, sdssSTD = tls.GetRV(sdssSampler)
    
    print "\n"
    print "sdss:"
    print sdssRV,sdssSTD
    print "\n"
    plot_format()
    for i in range(len(modelVels)):
        #plt.plot(sdssVels[i],sdssFluxes[i]+i*off)
        if i == 0:
            co = 'b'
        elif i == 1:
            co = 'g'
        else:
            co = 'r'
        plt.step(sdssVels[i],sdssFluxes[i]+i*off,where='mid',linewidth=1.5,color=co)
        
    plt.axvline(sdssRV,color='purple',ls='--')
    plt.axvline(sdssRV+sdssSTD,color='green',ls='--')
    plt.axvline(sdssRV-sdssSTD,color='green',ls='--')
    plt.axvline(0,color='k',ls='--')
    
    plt.plot(sdssVels[0],voigt(sdssVels[0],ld1,lw1,gd1,gw1,sdssRV)+0*off,color='k',linewidth=1.5,label='data fits')
    plt.plot(sdssVels[0],voigt(sdssVels[0],ld1,lw1,gd1,gw1,sdssRV+sdssSTD)+0*off,color='purple',linewidth=1.5,label='data fits + std',ls='--')
    plt.plot(sdssVels[0],voigt(sdssVels[0],ld1,lw1,gd1,gw1,sdssRV-sdssSTD)+0*off,color='purple',linewidth=1.5,label='data fits - std',ls='--')
    
    plt.plot(sdssVels[1],voigt(sdssVels[1],ld2,lw2,gd2,gw2,sdssRV)+1*off,color='k',linewidth=1.5)
    plt.plot(sdssVels[1],voigt(sdssVels[1],ld2,lw2,gd2,gw2,sdssRV+sdssSTD)+1*off,color='purple',linewidth=1.5,ls='--')
    plt.plot(sdssVels[1],voigt(sdssVels[1],ld2,lw2,gd2,gw2,sdssRV-sdssSTD)+1*off,color='purple',linewidth=1.5,ls='--')
    

    
    plt.plot(sdssVels[2],voigt(sdssVels[2],ld3,lw3,gd3,gw3,sdssRV)+2*off,color='k',linewidth=1.5)
    plt.plot(sdssVels[2],voigt(sdssVels[2],ld3,lw3,gd3,gw3,sdssRV+sdssSTD)+2*off,color='purple',linewidth=1.5,ls='--')
    plt.plot(sdssVels[2],voigt(sdssVels[2],ld3,lw3,gd3,gw3,sdssRV-sdssSTD)+2*off,color='purple',linewidth=1.5,ls='--')

    
    plt.xlim(-1500,1500)
    plt.title(wdName+" RV value="+str(sdssRV)+" RV Err="+str(sdssSTD))
    plt.savefig("sdssFits/"+wdName+"_sdssVelFit.pdf")
    plot_format()

    sdssVelFig = corner.corner(sdssSamplesChain,label=["RV"])
    sdssVelFig.savefig("sdssFits/corner_sdssVelFits.pdf")
    
    sampler = tls.MCMCfit(lnprobSum,args=(vels,fluxes,ferrs),nwalkers=nwalkers,ndim=ndim,burnInSteps=8000,steps=8000)


    ### Get the RV Shifts from the sampler with errors
    rvFit,rvStd = tls.GetRV(sampler)
    print rvFit,rvStd

    numArr.append(j)
    rvArr.append(rvFit)
    stdArr.append(rvStd)

    off = 0.5
    plot_format()
    ## VERTICAL LINES
    plt.axvline(0,ls="--",color='k')
    plt.axvline(rvFit,color='purple',linewidth=1.5,label="RV fit")
    plt.axvline(rvFit+rvStd,ls='--',color='purple',label="RV +/- $\sigma$")
    plt.axvline(rvFit-rvStd,ls='--',color='purple')

    ## DATA (MODEL AND REAL)
    for i in range(len(modelVels)):
        plt.plot(modelVels[i],modelFluxes[i]+i*off)
        if i == 0:
            co = 'b'
        elif i == 1:
            co = 'g'
        else:
            co = 'r'
        plt.step(vels[i],fluxes[i]+i*off,where='mid',linewidth=1.5,color=co)

    ### REAL DATA FITS
    plt.plot(vels[0],voigt(vels[0],ld1,lw1,gd1,gw1,rvFit)+0*off,color='k',linewidth=1.5,label='data fits')
    plt.plot(vels[1],voigt(vels[1],ld2,lw2,gd2,gw2,rvFit)+1*off,color='k',linewidth=1.5)
    plt.plot(vels[2],voigt(vels[2],ld3,lw3,gd3,gw3,rvFit)+2*off,color='k',linewidth=1.5)

    ### MODEL DATA FITS
    plt.plot(modelVels[0],voigtModel(modelVels[0],ld1,lw1,gd1,gw1)+0*off,color="y",linewidth=1.5,label='model fits')
    plt.plot(modelVels[1],voigtModel(modelVels[1],ld2,lw2,gd2,gw2)+1*off,color="y",linewidth=1.5)
    plt.plot(modelVels[2],voigtModel(modelVels[2],ld3,lw3,gd3,gw3)+2*off,color="y",linewidth=1.5)

    plt.xlim(min(vels[0]),max(vels[0]))
    plt.title(wdName+" Spectrum: " +str(j))
    plt.xlabel("velocity [km/s]")
    plt.ylabel("Normalized Flux + offset")
    plt.legend(prop={'size':8})
    plt.savefig("ModelFits/VelFits/"+wdName+"spec_"+str(j)+".pdf")
    plot_format()
                           
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
#plt.title(wdName + " velocity vs time\n{v[t]="+str(params[0])+"*sin[2*pi*(t/"+str(params[1])+") + "+str(params[2])+"] + "+str(params[3])+"}")
plt.savefig("ModelFits/"+wdName+"_time.pdf")

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
plt.title(wdName + " velocity vs time"+'\n'+"$v(t) = {0:.3f}*sin(2\pi(t/{1:.3f}) + {2:.3f}) + {3:.3f}$".format(params[0],params[1],params[2],params[3]))

plt.subplot(4,1,2)
plt.errorbar(timeArr,rvArr,yerr=stdArr,linestyle='None',marker='o')
#for A,P,Ph,Gam in samplesChain[np.random.randint(len(samplesChain),size=100)]:
    #plt.plot(newTime,sine(newTime,A,P,Ph,Gam),color='k',alpha=0.1)
    #print A,P,Ph,Gam
plt.plot(newTime,sine(newTime,params[0],params[1],params[2],params[3]),color="red",alpha=0.75)
plt.xlabel("MJD [days]")
plt.ylabel("Velocity [km/s]")
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
#plt.title(wdName + " velocity vs time")
plt.xlim(np.max(timeArr)-0.1,np.max(timeArr)+0.1)

plt.savefig("ModelFits/"+wdName+"_time_zoomed.pdf")

#print ""
#print Amp,Per,Phi,Gamma

fig = corner.corner(samplesChain, labels=["A","P","Phi","Gam"])
fig.savefig("ModelFits/"+wdName+"_Triangle.pdf")


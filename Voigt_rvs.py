import tools as tls
from astropy.io import fits
import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt
from plot_format import plot_format
import os
import sys
import emcee as mc


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

lines = [line.rstrip('\n') for line in open('filelist')]
plot_format()
#print lines
plotNumber = 1
plotIndex = 1

for j in range(len(lines)):
#for j in range(0,1):

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
        ferr = errorNorm[wherr]
        wl = np.linspace(lowerLine,upperLine,len(flux))
        #plt.plot(wl,flux)
        #plt.show()
    
        vel = []
        for w in range(len(wl)):
            v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
            vel.append(v)
        if (lineIndex == 1):
            velBeta = vel
            fluxBeta = flux
            ferrBeta = ferr
        if (lineIndex == 2):
            velGamma = vel
            fluxGamma = flux
            ferrGamma = ferr
        if (lineIndex == 3):
            velDelta = vel
            fluxDelta = flux
            ferrDelta = ferr
    ############################
    ###### FITTING STUFF #######
    ############################
        
    vels = [velBeta,velGamma,velDelta]
    fluxes = [fluxBeta,fluxGamma,fluxDelta]
    ferrs = [ferrBeta,ferrGamma,ferrDelta]

    ndim, nwalkers = 13, 100
    
    sampler = mc.EnsembleSampler(nwalkers,ndim,lnprobSum,args=(vels,fluxes,ferrs))
    np.random.seed(1234)
    p0  = np.random.random((nwalkers, ndim))
    
    
    print "running burn in"
    p0, prob, _ = sampler.run_mcmc(p0,4000)
    sampler.reset()
    
    print "running mcmc production"
    p0, _, _ = sampler.run_mcmc(p0,4000)
    flatchain = sampler.flatchain
    samples = sampler.flatchain.reshape(-1,ndim)
    samplesTrans = sampler.flatchain.reshape(-1,ndim).T
    
    ld1FitArr = samplesTrans[0]
    lw1FitArr = samplesTrans[1]
    gd1FitArr = samplesTrans[2]
    gw1FitArr = samplesTrans[3]
    ld2FitArr = samplesTrans[4]
    lw2FitArr = samplesTrans[5]
    gd2FitArr = samplesTrans[6]
    gw2FitArr = samplesTrans[7]
    ld3FitArr = samplesTrans[8]
    lw3FitArr = samplesTrans[9]
    gd3FitArr = samplesTrans[10]
    gw3FitArr = samplesTrans[11]
    RVFitArr = samplesTrans[12]

    ld1Fit = ld1FitArr.mean()
    lw1Fit = lw1FitArr.mean()
    gd1Fit = gd1FitArr.mean()
    gw1Fit = gw1FitArr.mean()
    ld2Fit = ld2FitArr.mean()
    lw2Fit = lw2FitArr.mean()
    gd2Fit = gd2FitArr.mean()
    gw2Fit = gw2FitArr.mean()
    ld3Fit = ld3FitArr.mean()
    lw3Fit = lw3FitArr.mean()
    gd3Fit = gd3FitArr.mean()
    gw3Fit = gw3FitArr.mean()
    RVFit = RVFitArr.mean()

    ld1Std = ld1FitArr.std()
    lw1Std = lw1FitArr.std()
    gd1Std = gd1FitArr.std()
    gw1Std = gw1FitArr.std()
    ld2Std = ld2FitArr.std()
    lw2Std = lw2FitArr.std()
    gd2Std = gd2FitArr.std()
    gw2Std = gw2FitArr.std()
    ld3Std = ld3FitArr.std()
    lw3Std = lw3FitArr.std()
    gd3Std = gd3FitArr.std()
    gw3Std = gw3FitArr.std()
    RVStd = RVFitArr.std()


    #Plot Hbeta vertical line
    #plt.figure(1).add_axes((.1,.3,.8,.6))
    plt.subplot(2,2,plotIndex)
    plt.plot(Owl,Normflux)
    plt.axhline(1,ls='--',color='black')
    plt.subplot(2,2,plotIndex+1)
    plt.axvline(0,ls='--',color='black')
    plt.axvline(RVFit,ls='--',color='purple',linewidth=1.5)
    off = 0.3
    
    plt.step(vels[0],fluxes[0]+0*off,where='mid',linewidth=2,label=lineNames[1]+" - "+str(RVStd))
    plt.step(vels[1],fluxes[1]+1*off,where='mid',linewidth=2,label=lineNames[2]+" - "+str(RVStd))
    plt.step(vels[2],fluxes[2]+2*off,where='mid',linewidth=2,label=lineNames[3]+" - "+str(RVStd))

    plt.plot(vels[0],voigt(vels[0],ld1Fit,lw1Fit,gd1Fit,gw1Fit,RVFit)+0*off,color='black',linewidth=2)
    plt.plot(vels[1],voigt(vels[1],ld2Fit,lw2Fit,gd2Fit,gw2Fit,RVFit)+1*off,color='black',linewidth=2)
    plt.plot(vels[2],voigt(vels[2],ld3Fit,lw3Fit,gd3Fit,gw3Fit,RVFit)+2*off,color='black',linewidth=2)

    plt.xlim(min(vel),max(vel))
    plt.title("Spectrum: "+str(j))
    plt.xlabel("velocity [km/s]")
    plt.ylabel("Normalized flux + offset")
    plt.legend(prop={'size':8})
    
    
    plotIndex += 2    
    if plotIndex == 5:
        plt.savefig("VoigtFitPlots/"+wdName+"_fit_pg"+str(plotNumber)+".pdf")
        plot_format()
        plotNumber += 1
        plotIndex = 1
    

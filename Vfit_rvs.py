import tools as tls
from astropy.io import fits
import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt
from plot_format import plot_format
import os
import sys
import emcee as mc
import corner as cn

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

    
    ndim, nwalkers = 7, 100
    sampler = mc.EnsembleSampler(nwalkers,ndim,lnprobSum,args=(vels,fluxes,ferrs))
    np.random.seed(1234)
    p0  = np.random.random((nwalkers, ndim))

    p0, prob, _ = sampler.run_mcmc(p0,500)
    sampler.reset()

    p0, _, _ = sampler.run_mcmc(p0,1000)
    flatchain = sampler.flatchain
    samples = sampler.flatchain.reshape(-1,ndim).T

    d1FitArr = samples[0]
    w1FitArr = samples[1]
    d2FitArr = samples[2]
    w2FitArr = samples[3]
    d3FitArr = samples[4]
    w3FitArr = samples[5]
    RVFitArr = samples[6]
    
    d1Fit = d1FitArr.mean()
    w1Fit = w1FitArr.mean()
    d2Fit = d2FitArr.mean()
    w2Fit = w2FitArr.mean()
    d3Fit = d3FitArr.mean()
    w3Fit = w3FitArr.mean()


    RVFit = RVFitArr.mean()
    RVStd = RVFitArr.std()
    off = 0.3

    #if (plotIndex == 1):
    #    plt.hist(RVFitArr)
    #    plt.title("RV values (Fit: "+str(RVFit)+") (std: "+str(RVStd)+")")
    #    plt.savefig("cornerPlots/"+wdName+"_corner_"+str(j)+".pdf")
    #    plot_format()

    plt.subplot(2,2,plotIndex)
    plt.plot(Owl,Normflux)
    plt.axhline(1,ls='--',color='black')
    plt.subplot(2,2,plotIndex+1)
    plt.axvline(0,ls='--',color='black')
    plt.axvline(RVFit,ls='--',color='purple',linewidth=1.5)

    plt.step(vels[0],fluxes[0]+0*off,where='mid',linewidth=2,label=lineNames[1]+" - "+str(RVStd))
    plt.step(vels[1],fluxes[1]+1*off,where='mid',linewidth=2,label=lineNames[2]+" - "+str(RVStd))
    plt.step(vels[2],fluxes[2]+2*off,where='mid',linewidth=2,label=lineNames[3]+" - "+str(RVStd))
    
    plt.plot(vels[0],lorentzian(vels[0],d1Fit,w1Fit,RVFit)+0*off,color='black',linewidth=2)
    plt.plot(vels[1],lorentzian(vels[1],d2Fit,w2Fit,RVFit)+1*off,color='black',linewidth=2)
    plt.plot(vels[2],lorentzian(vels[2],d3Fit,w3Fit,RVFit)+2*off,color='black',linewidth=2)

    plt.xlim(min(vel),max(vel))
    plt.title("Spectrum: "+str(j))
    plt.xlabel("velocity [km/s]")
    plt.ylabel("Normalized flux + offset")
    plt.legend(prop={'size':5})

    plotIndex += 2    
    if plotIndex == 5:
        plt.savefig("fitPlots/"+wdName+"_fit_pg"+str(plotNumber)+".pdf")
        plot_format()
        plotNumber += 1
        plotIndex = 1
        plt.hist(RVFitArr)
        plt.title("RV values (Fit: "+str(RVFit)+") (std: "+str(RVStd)+")")
        plt.savefig("cornerPlots/"+wdName+"_corner_"+str(j)+".pdf")
        plot_format()

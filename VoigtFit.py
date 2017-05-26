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
    #print header
    print header['UTMIDDLE']
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

        ############################
        ###### FITTING STUFF #######
        ############################
    
        ndim, nwalkers = 5, 100
    
        sampler = mc.EnsembleSampler(nwalkers,ndim,lnprob,args=(vel,flux,ferr))
        np.random.seed(1234)
        p0  = np.random.random((nwalkers, ndim))
        
        
        print "running burn in"
        p0, prob, _ = sampler.run_mcmc(p0,1000)
        sampler.reset()
        
        print "running mcmc production"
        p0, _, _ = sampler.run_mcmc(p0,2000)
        flatchain = sampler.flatchain
        samples = sampler.flatchain.reshape(-1,ndim)
        samplesTrans = sampler.flatchain.reshape(-1,ndim).T

        LdepthFitArr = samplesTrans[0]
        LwidthFitArr = samplesTrans[1]
        GdepthFitArr = samplesTrans[2]
        GwidthFitArr = samplesTrans[3]
        RVFitArr = samplesTrans[4]
        
        LdepthMean = LdepthFitArr.mean()
        LwidthMean = LwidthFitArr.mean()
        GdepthMean = GdepthFitArr.mean()
        GwidthMean = GwidthFitArr.mean()
        RVMean = RVFitArr.mean()
    
        LdepthStd = LdepthFitArr.std()
        LwidthStd = LwidthFitArr.std()
        GdepthStd = GdepthFitArr.std()
        GwidthStd = GwidthFitArr.std()
        RVStd = RVFitArr.std()

        LdepthFit = LdepthMean
        LwidthFit = LwidthMean
        GdepthFit = GdepthMean
        GwidthFit = GwidthMean
        RVFit = RVMean

        #Plot Hbeta vertical line
    #plt.figure(1).add_axes((.1,.3,.8,.6))
        plt.subplot(2,2,plotIndex)
        plt.plot(Owl,Normflux)
        plt.axhline(1,ls='--',color='black')
        plt.subplot(2,2,plotIndex+1)
        plt.axvline(0,ls='--',color='black')
        plt.axvline(RVFit,ls='--',color='purple',linewidth=1.5)
        off = 0.3
    
        plt.step(vel,flux+lineIndex*off,where='mid',linewidth=2,label=lineNames[lineIndex] +" - " + str(RVStd))
        plt.plot(vel,voigt(vel,LdepthFit,LwidthFit,GdepthFit,GwidthFit,RVFit)+lineIndex*off,color='black',linewidth=2)
    
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
    


"""    import corner as cn

    fig = cn.corner(flatchain,labels=["depth","width","RV shift"])
    fig.savefig("fitPlots/"+wdName+"_corner.pdf")
"""

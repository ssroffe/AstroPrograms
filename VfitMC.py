import tools as tls
from astropy.io import fits
import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt
from plot_format import plot_format
import os
import sys
import emcee as mc


def lorentzian(x, depth, width, RVShift):
    return 1.0-depth/(1.0 + ((x-RVShift)/width)**2)

def lnlike(p,x,y,err):
    depth,width,RVShift = p
    return -np.sum((y-lorentzian(x,depth,width,RVShift))**2/(2*err**2))

def lnprior(p):
    depth,width,RVShift = p
    if 0.2 < depth < 3.0 and 0.0 < width < 1000 and -1000.0 < RVShift < 1000.0:
        return 0.0
    return -np.inf

def lnprob(p, x, y, yerr):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p, x, y, yerr)

lines = [line.rstrip('\n') for line in open('filelist')]
#plot_format()
#print lines
plotNumber = 1
plotIndex = 1

#for j in range(len(lines)):
for j in range(0,1):

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
    #for lineIndex in range(1,4):
    for lineIndex in range(1,2):
        
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
    
        ndim, nwalkers = 3, 100
    
        sampler = mc.EnsembleSampler(nwalkers,ndim,lnprob,args=(vel,flux,ferr))
        np.random.seed(1234)
        p0  = np.random.random((nwalkers, ndim))
        
        #p0[0] = p0[0] + 0.1
        #p0 = [[0.1,0.1,0.1] for i in range(nwalkers)]
        
        print "running burn in"
        p0, prob, _ = sampler.run_mcmc(p0,1000)
        sampler.reset()
        
        print "running mcmc production"
        p0, _, _ = sampler.run_mcmc(p0,2000)
        flatchain = sampler.flatchain
        samples = sampler.flatchain.reshape(-1,ndim).T

        depthFitArr = samples[0]
        widthFitArr = samples[1]
        RVFitArr = samples[2]
        
        depthMean = depthFitArr.mean()
        widthMean = widthFitArr.mean()
        RVMean = RVFitArr.mean()
    
        depthStd = depthFitArr.std()
        widthStd = widthFitArr.std()
        RVStd = RVFitArr.std()

        depthFit = depthMean
        widthFit = widthMean
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
    
        plt.step(vel,flux+lineIndex*off,where='mid',linewidth=2,label=lineNames[lineIndex])
        plt.plot(vel,lorentzian(vel,depthFit,widthFit,RVFit)+lineIndex*off,color='black',linewidth=2)
    
        plt.xlim(min(vel),max(vel))
        plt.title("Spectrum: "+str(j))
        plt.xlabel("velocity [km/s]")
        plt.ylabel("Normalized flux + offset")
        plt.legend(prop={'size':5})


    import corner as cn

    fig = cn.corner(flatchain,labels=["depth","width","RV shift"])
    fig.savefig("fitPlots/"+wdName+"_corner.pdf")

"""    
    plotIndex += 2    
    if plotIndex == 5:
        plt.savefig("fitPlots/"+wdName+"_fit_pg"+str(plotNumber)+".pdf")
        plot_format()
        plotNumber += 1
        plotIndex = 1
    
"""

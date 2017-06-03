#################################
##### ORBITAL FITTING TOOLS #####
#################################

### Make a directory
def mkdir(directory):
    import os

    if not os.path.exists(directory):
        os.makedirs(directory)
    return

### Get filename
def GetFileName(path):
    import os
    return os.path.basename(path)[:-5]

### Get the middle of the observation in UTC
def GetDateTime(path):
    from astropy.io import fits
    import numpy as np

    header = fits.getheader(path)
    #print header
    dateTime = header['UTMIDDLE']
    return dateTime

### Get JD of the start of the observation
def GetJD(path):
    from astropy.io import fits
    
    header = fits.getheader(path)
    #print header
    jd = header['JD']
    return jd

### Formatting plots to look nice
def PlotFormat():
    #Plots
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.rcParams.update({'figure.autolayout':'True'})
    plt.rcParams.update({'font.size': 12})
    plt.rcParams.update({'mathtext.default':'regular'})
    plt.rcParams.update({'mathtext.fontset':'stixsans'})
    plt.rcParams.update({'axes.linewidth': 1.5})
    plt.rcParams.update({'xtick.major.size': 5})
    plt.rcParams.update({'xtick.major.width': 1.25 })
    plt.rcParams.update({'xtick.minor.size': 2.5})
    plt.rcParams.update({'xtick.minor.width': 1.25 })
    plt.rcParams.update({'ytick.major.size': 5})
    plt.rcParams.update({'ytick.major.width': 1.25 })
    plt.rcParams.update({'ytick.minor.size': 2.5})
    plt.rcParams.update({'ytick.minor.width': 1.25 })
    plt.rc('legend',**{'fontsize':'x-small'})
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()

    return 0

### Get the raw spectrum
def RawSpectrum(path):
    from astropy.io import fits
    import numpy as np
    
    header = fits.getheader(path)

    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    flux = fits.getdata(path)
    wl = np.linspace(wl_start,wl_end,len(flux))

    return (wl,flux)


### Normalize the spectrum to 1 and plot it
def Normalize(path):
    from astropy.io import fits
    import numpy as np
    import scipy.interpolate as sp
    import matplotlib.pyplot as plt

    basename = GetFileName(path)
    
    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    flux = np.array(fits.getdata(path))
    wl = np.linspace(wl_start,wl_end,len(flux))

    Oflux = np.array(fits.getdata(path))
    Owl = np.linspace(wl_start,wl_end,len(Oflux))

    flux = flux / max(flux)
    Oflux = Oflux / max(Oflux)
    
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])

    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])

    upperLine = lineList + 50
    lowerLine = lineList - 50

    ## cuts
    for i in range(len(lineList)):
        cutList = np.where((wl >= lowerLine[i]) & (wl <= upperLine[i]))
        flux = np.delete(flux,cutList)
        wl = np.delete(wl,cutList)
    
    ## interpolation
    interp = sp.splrep(wl,flux,s=5,k=3)

    wlnew = np.linspace(wl_start,wl_end,len(Oflux))
    fluxnew = sp.splev(wlnew,interp,der=0)

    normalization = Oflux/fluxnew

    errorPath = path[:-5] + "variance.fits"

    errorFlux = fits.getdata(errorPath)
    errorWl = np.linspace(wl_start,wl_end,len(errorFlux))
    errorNorm = errorFlux / fluxnew

    ## Plot things

    # Error plots
    
    PlotFormat()
    plt.plot(errorWl,errorFlux,'b',alpha=0.5,label="Original")
    plt.plot(wlnew,fluxnew,'g',alpha=0.7,label="Interpolated")
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Flux")
    plt.title(basename+" Error")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/RawErrorPlots/"+basename+"_error.pdf")
    plt.figure()

    PlotFormat()
    plt.plot(errorWl,errorNorm,'r',alpha=0.3,label="Normalized")
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Flux")
    plt.title(basename+" Error [Normalized]")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/ErrorPlots/"+basename+"_normerror.pdf")
    plt.figure()

    # Raw stuff
    PlotFormat()
    plt.plot(Owl,Oflux,'b',alpha=0.25,label="original")
    plt.plot(wl,flux,'r',alpha=0.3,label="cut")
    plt.plot(wlnew,fluxnew,'g',alpha=0.7,label="interpolated")
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux")
    plt.title(basename)
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/RawPlots/"+basename+".pdf")
    plt.figure()

    # Normalized
    PlotFormat()
    plt.plot(Owl,normalization,'b',ls='-',alpha=0.4,label="normalized")
    plt.plot(Owl,np.ones(len(Owl)),'black',ls='--',alpha=0.75)
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux [Normalized]")
    plt.title(basename+" [Normalized]")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/plots/"+basename+"_Norm.pdf")
    plt.figure()

    # Combine Original and normalized
    PlotFormat()
    plt.subplot(2,1,1)
    plt.plot(Owl,normalization,'b',ls='-',alpha=0.4,label="normalized")
    plt.plot(Owl,np.ones(len(Owl)),'black',ls='--',alpha=0.75)
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux [Normalized]")
    plt.ylim(-1.0,3.0)
    plt.title(basename+" [Normalized]")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    
    plt.subplot(2,1,2)
    plt.plot(Owl,Oflux,'b',alpha=0.25,label="original")
    plt.plot(wl,flux,'r',alpha=0.3,label="cut")
    plt.plot(wlnew,fluxnew,'g',alpha=0.7,label="interpolated")
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux")
    plt.title(basename)
    plt.xlim(wl_start,wl_end)
    plt.legend()

    plt.savefig("/home/seth/research/CombinedPlots/"+basename+"_combined.pdf")
    
    return (Owl,normalization,errorNorm)


### Normalize the spectrum to 1 and plot it with SG filter
def NormalizeSG(path):
    from astropy.io import fits
    import numpy as np
    import scipy.interpolate as sp
    import matplotlib.pyplot as plt

    basename = GetFileName(path)
    
    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    flux = fits.getdata(path)
    wl = np.linspace(wl_start,wl_end,len(flux))

    Oflux = fits.getdata(path)
    Owl = np.linspace(wl_start,wl_end,len(Oflux))

    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])

    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    
    upperLine = lineList + 50
    lowerLine = lineList - 50

    ## cuts
    for i in range(len(lineList)):
        cutList = np.where((wl >= lowerLine[i]) & (wl <= upperLine[i]))
        flux = np.delete(flux,cutList)
        wl = np.delete(wl,cutList)
    
    ## interpolation
    interp = sp.splrep(wl,flux,s=5,k=3)
    
    wlnew = np.linspace(wl_start,wl_end,len(Oflux))
    fluxnew = sp.splev(wlnew,interp,der=0)

    #rft = np.fft.rfft(fluxnew)
    #rft[5:] = 0
    #fluxnew = np.fft.irfft(rft)

    #fluxnew = savitzky_golay(fluxnew,51,3)

    #print fluxnew
    #plt.plot(wlnew,fluxnew)
    #plt.show()
    #plt.clf()

    normalization = Oflux/fluxnew
    normalizationNew = savitzky_golay(normalization,103,5)
    normalizationNew = ((normalizationNew - 1)*2)+1
    #rft = np.fft.rfft(normalization)
    #normalization = np.fft.irfft(rft)

    errorPath = path[:-5] + "variance.fits"

    errorFlux = fits.getdata(errorPath)
    errorWl = np.linspace(wl_start,wl_end,len(errorFlux))
    errorNorm = errorFlux / fluxnew

    ## Plot things

    # Error plots
    
    PlotFormat()
    plt.plot(errorWl,errorFlux,'b',alpha=0.5,label="Original")
    plt.plot(wlnew,fluxnew,'g',alpha=0.7,label="Interpolated")
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Flux")
    plt.title(basename+" Error")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/RawErrorPlots/"+basename+"_error.pdf")
    plt.figure()

    PlotFormat()
    plt.plot(errorWl,errorNorm,'r',alpha=0.3,label="Normalized")
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Flux")
    plt.title(basename+" Error [Normalized]")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/ErrorPlots/"+basename+"_normerror.pdf")
    plt.figure()

    # Raw stuff
    PlotFormat()
    plt.plot(Owl,Oflux,'b',alpha=0.25,label="original")
    plt.plot(wl,flux,'r',alpha=0.3,label="cut")
    plt.plot(wlnew,fluxnew,'g',alpha=0.7,label="interpolated")
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux")
    plt.title(basename)
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/RawPlots/"+basename+".pdf")
    plt.figure()

    # Normalized
    PlotFormat()
    plt.plot(Owl,normalization,'b',ls='-',alpha=0.4,label="normalized")
    plt.plot(Owl,normalizationNew,'r',ls='-',alpha=0.4,label="normalized")
    for i in lineList:
        plt.axvline(i,ls='--',color='k')
    plt.plot(Owl,np.ones(len(Owl)),'black',ls='--',alpha=0.75)
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux [Normalized]")
    plt.title(basename+" [Normalized]")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    plt.savefig("/home/seth/research/plots/"+basename+"_Norm.pdf")
    plt.figure()

    # Combine Original and normalized
    PlotFormat()
    plt.subplot(2,1,1)
    plt.plot(Owl,normalization,'b',ls='-',alpha=0.4,label="normalized")
    plt.plot(Owl,np.ones(len(Owl)),'black',ls='--',alpha=0.75)
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux [Normalized]")
    plt.ylim(-1.0,3.0)
    plt.title(basename+" [Normalized]")
    plt.xlim(wl_start,wl_end)
    plt.legend()
    
    plt.subplot(2,1,2)
    plt.plot(Owl,Oflux,'b',alpha=0.25,label="original")
    plt.plot(wl,flux,'r',alpha=0.3,label="cut")
    plt.plot(wlnew,fluxnew,'g',alpha=0.7,label="interpolated")
    plt.xlabel("wavelength [$\AA$]")
    plt.ylabel("flux")
    plt.title(basename)
    plt.xlim(wl_start,wl_end)
    plt.legend()

    plt.savefig("/home/seth/research/CombinedPlots/"+basename+"_combined.pdf")
    
    return (Owl,normalization,errorNorm)

## Normalize spectrum without Plotting
def NormNoPlot(path):
    from astropy.io import fits
    import numpy as np
    import scipy.interpolate as sp
    import matplotlib.pyplot as plt

    basename = GetFileName(path)
    
    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))

    flux = np.array(fits.getdata(path))
    wl = np.linspace(wl_start,wl_end,len(flux))
    
    Oflux = np.array(fits.getdata(path))
    Owl = np.linspace(wl_start,wl_end,len(Oflux))

    #flux = flux / max(flux)
    #Oflux = Oflux / max(Oflux)
    
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[3975.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])

    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4030.0, 4180.0], [4250.0,4475.0], [4700.0,4980.0]])
    
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    offset = 100
    upperLine = lineList + offset
    lowerLine = lineList - offset

    # cuts
    for i in range(len(lineList)):
        cutList = np.where((wl >= lineWindows[i][0]) & (wl <= lineWindows[i][1]))
        #cutList = np.where((wl >= lowerLine[i]) & (wl <= upperLine[i]))
        flux = np.delete(flux,cutList)
        wl = np.delete(wl,cutList)
    
    ## interpolation
    #interp = sp.splrep(wl,flux,s=5,k=3)

    #wlnew = np.linspace(wl_start,wl_end,len(Oflux))
    #fluxnew = sp.splev(Owl,interp,der=0)
    #fluxnew = savitzky_golay(flux,53,5)

    interp = np.interp(Owl,wl,flux)
    
    fluxnew = savitzky_golay(interp,103,5)
    
    #plt.plot(Owl,fluxnew,alpha=0.7)
    #plt.plot(Owl,interp)
    #plt.plot(Owl,Oflux,alpha=0.4)
    #plt.title("TMP2")
    #plt.show()
    normalization = Oflux/fluxnew

    ####NEW STUFF

    errorPath = path[:-5] + "variance.fits"

    errorFlux = fits.getdata(errorPath)
    errorWl = np.linspace(wl_start,wl_end,len(errorFlux))
    errorNorm = errorFlux / fluxnew

    #normalization = ((normalization - 1)*2)+1
    #errorNorm = ((errorNorm - 1)*2)+1
    
    return (Owl,normalization,errorNorm)

#Function to plot chains
def plotchains(sampler,nParamsPerPage,filePrefix,fileSuffix):

    import matplotlib.pyplot as plt
    import numpy as np
    nParams =  (sampler.chain.shape)[2]
    nWalkers =  (sampler.chain.shape)[0]
   
    plt.clf()
    plt.figure(figsize=(20, 20))
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    nPages = np.ceil(nParams/nParamsPerPage)
    thisPage = 0

    for parameter in range(nParams):

        parName = '$\Psi_{' + str(parameter) + '}$'
   
        indexPlot = 2*(parameter-thisPage*nParamsPerPage)+1
        plt.subplot(nParamsPerPage, 2, indexPlot)
        plt.xlabel('Step')
        plt.ylabel(parName)
        for i in range(nWalkers): plt.plot(sampler.chain[i,:,parameter])
       
        indexPlot += 1
        plt.subplot(nParamsPerPage, 2, indexPlot)
        plt.xlabel(parName)
        plt.ylabel('N')
        plt.hist(sampler.flatchain[:,parameter],bins = 20)

        if (np.mod(parameter,nParamsPerPage) == nParamsPerPage-1) or (parameter == nParams-1) :
            pageTitle = '_Params_'+str(parameter-nParamsPerPage+1)+'_'+str(parameter)
            plt.savefig(filePrefix+pageTitle+fileSuffix)
            thisPage = thisPage + 1
            plt.clf()
            plt.figure(figsize=(20, 20))

### Get the velocities of a single H line
def GetVelArr(path,lineIndex):
    import numpy as np
    from astropy.io import fits
    import os

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    timeTaken = basename[15:]
    #timeTaken = basename[7:]
    #print timeTaken

    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    
    Owl,Normflux,errorNorm = NormNoPlot(path)
    
    ## Amplify the lines
    #Normflux = ((Normflux-1)*2)+1

    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineWindows = np.array([[4000.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])

    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"])

    print("Obtaining "+lineNames[lineIndex]+" velocities")

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
    wl = Owl[wherr]
    #plt.plot(wl,flux)
    #plt.show()
    
    vel = []
    for w in range(len(wl)):
        v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
        vel.append(v)

    return (vel,flux,ferr)

### Get the velocities of all three H lines (beta,gamma,delta)
def GetAllVelocities(path):
    import numpy as np
    from astropy.io import fits
    import os

    c = 299792.458 #km/s
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    timeTaken = basename[15:]
    #timeTaken = basename[7:]
    #print timeTaken

    header = fits.getheader(path)
    wl_start, wl_end = map(float,header['W_RANGE'].split(" "))
    
    Owl,Normflux,errorNorm = NormNoPlot(path)
    
    ## Amplify the lines
    #Normflux = ((Normflux-1)*2)+1

    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4000.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])

    lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    lineWindows = np.array([[4030.0, 4180.0], [4250.0,4475.0], [4700.0,4980.0]])
    
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"])
    #lineIndex = 1
    #for lineIndex in range(1,4):
    for lineIndex in range(len(lineList)):
        
        #offset = 50
        #offset = 25
        offset = 30
    
        #upperLine = lineWindows[lineIndex][1]
        #lowerLine = lineWindows[lineIndex][0]

        upperLine = lineList[lineIndex] + offset
        lowerLine = lineList[lineIndex] - offset
        
        #plt.axvline(upperLine,color='black')
        #plt.axvline(lowerLine,color="black")
        #plt.plot(Owl,Normflux)
        #plt.show()
    
        wherr = np.where((Owl >= lowerLine) & (Owl <= upperLine))
        flux = Normflux[wherr]
        ferr = errorNorm[wherr]
        wl = Owl[wherr]
        #wl = np.linspace(lowerLine,upperLine,len(flux))
        #plt.plot(wl,flux)
        #plt.show()
    
        vel = []
        for w in range(len(wl)):
            v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
            vel.append(v)
        if (lineIndex == 0):
            velBeta = vel
            fluxBeta = flux
            ferrBeta = ferr
        if (lineIndex == 1):
            velGamma = vel
            fluxGamma = flux
            ferrGamma = ferr
        if (lineIndex == 2):
            velDelta = vel
            fluxDelta = flux
            ferrDelta = ferr

    vels = np.array([np.array(velBeta),np.array(velGamma),np.array(velDelta)])
    fluxes = np.array([fluxBeta,fluxGamma,fluxDelta])
    ferrs = np.array([ferrBeta,ferrGamma,ferrDelta])

    return (vels,fluxes,ferrs)

### Do an MCMC fit
def MCMCfit(lnprob,args,burnInSampler=False,nwalkers=100,ndim=3,burnInSteps=1000,steps=2000,p=None):
    ## lnprob = the function of fitting with priors
    ## args = ntuple of the arguments e.g. (vel,flux,ferr)
    
    import emcee as mc
    import numpy as np
    import random as rand

    sampler = mc.EnsembleSampler(nwalkers,ndim,lnprob,args=args)
    np.random.seed(1234)

    if p == None:
        p0  = np.random.random((nwalkers, ndim)) / 10000 + 0.039
    else:
        p0 = p
    # burn in
    p0, _, _ = sampler.run_mcmc(p0, burnInSteps)

    if burnInSampler == True:
        return sampler

    sampler.reset()
    # production
    p0, _, _ = sampler.run_mcmc(p0,steps)

    return sampler

### Get the Radial Velocities from the MCMC fit
def GetRV(sampler):
    import numpy as np
    
    flatchain = sampler.flatchain
    ndim = np.shape(flatchain)[1]
    samplesChain = sampler.chain[:,:,:].reshape((-1,ndim))
    samples = sampler.flatchain.reshape(-1,ndim).T

    RVFitArr = samples[-1]
    rvMCMC = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samplesChain, [16,50,84],axis=0)))[-1]
    RVFit = RVFitArr.mean()
    RVStd = RVFitArr.std()

    #import corner
    #fig = corner.corner(samplesChain)
    #fig.show()
    #fig.savefig("/home/seth/testcorner.pdf")

    #print rvMCMC
    #print RVFit,RVStd
    
    return (RVFit, RVStd)

### Savitzky Golay Filter 
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def ModelNormNoPlot(path):
    from astropy.io import fits
    import numpy as np
    import scipy.interpolate as sp
    import matplotlib.pyplot as plt

    modelData = np.genfromtxt(path,skip_header=34)
    flux = np.array(modelData[:,1])
    wl = modelData[:,0]
    
    Oflux = np.array(modelData[:,1])
    Owl = modelData[:,0]

    #flux = flux / max(flux)
    #Oflux = Oflux / max(Oflux)
    
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4000.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    
    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4030.0, 4180.0], [4250.0,4475.0], [4700.0,4980.0]])
    
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    offset = 100
    upperLine = lineList + offset
    lowerLine = lineList - offset

    ## cuts
    for i in range(len(lineList)):
        cutList = np.where((wl >= lineWindows[i][0]) & (wl <= lineWindows[i][1]))
        #cutList = np.where((wl >= lowerLine[i]) & (wl <= upperLine[i]))
        flux = np.delete(flux,cutList)
        wl = np.delete(wl,cutList)

        
    ## interpolation
    #interp = sp.splrep(wl,flux,s=5,k=5)

    #fluxnew = sp.splev(Owl,interp,der=0)
    interp = np.interp(Owl,wl,flux)
    fluxnew = savitzky_golay(interp,53,5)


    #interp = sp.splrep(wl,fluxsg,s=0)
    #fluxnew = sp.interp1d(Owl,fluxsg)
    #fluxnew = sp.splev(Owl,interp,der=0)

    #interp = sp.InterpolatedUnivariateSpline(wl,flux,k=2)
    #fluxnew = interp(Owl)

    #fluxnew = np.interp(Owl,wl,flux)
    
    #plt.plot(Owl,Oflux)
    #plt.plot(Owl,fluxnew)
    #plt.title("TMP")
    #plt.show()
    
    normalization = Oflux/fluxnew

    ####NEW STUFF
    #normalization = ((normalization - 1)*2)+1
    
    return (Owl,normalization)

def ModelGetAllVelocities(path):
    import numpy as np
    from astropy.io import fits
    import os
    import matplotlib.pyplot as plt

    c = 299792.458 #km/s
    
    Owl,Normflux = ModelNormNoPlot(path)
    #plt.plot(Owl,Normflux)
    #plt.show()
    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    #lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4000.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])

    lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    lineWindows = np.array([[4030.0, 4180.0], [4250.0,4475.0], [4700.0,4980.0]])    
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"]
    )
    for lineIndex in range(len(lineList)):
        #offset = 50
        #offset = 25
        #offset = 30
        
        upperLine = lineWindows[lineIndex][1]
        lowerLine = lineWindows[lineIndex][0]

        #upperLine = lineList[lineIndex] + offset
        #lowerLine = lineList[lineIndex] - offset
        
        #plt.axvline(upperLine,color='black')
        #plt.axvline(lowerLine,color="black")
        #plt.plot(Owl,Normflux)
        #plt.show()
    
        wherr = np.where((Owl >= lowerLine) & (Owl <= upperLine))
        flux = Normflux[wherr]
        #ferr = errorNorm[wherr]
        wl = Owl[wherr]
        #wl = np.linspace(lowerLine,upperLine,len(flux))
        #plt.plot(wl,flux)
        #plt.show()
    
        vel = []
        for w in range(len(wl)):
            v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
            vel.append(v)
        if (lineIndex == 0):
            velBeta = vel
            fluxBeta = flux
            #ferrBeta = ferr
        if (lineIndex == 1):
            velGamma = vel
            fluxGamma = flux
            #ferrGamma = ferr
        if (lineIndex == 2):
            velDelta = vel
            fluxDelta = flux
            #ferrDelta = ferr

    vels = np.array([np.array(velBeta),np.array(velGamma),np.array(velDelta)])
    fluxes = np.array([fluxBeta,fluxGamma,fluxDelta])

    return (vels,fluxes)

def CSVNormNoPlot(path):
    from astropy.io import fits
    import numpy as np
    import scipy.interpolate as sp
    import matplotlib.pyplot as plt

    sdssData = np.genfromtxt(path,skip_header=1,delimiter=',')
    flux = np.array(sdssData[:,1])
    wl = np.array(sdssData[:,0])

    Teff = 11980
    
    Oflux = np.array(sdssData[:,1])

    #flux = flux /max(flux)
    #Oflux = Oflux / max(Oflux)
    
    Owl = sdssData[:,0]
    #offset = 75
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    #lineWindows = np.array([[6563.0-offset,6563.0+offset],[4800.0,4950.0],[4300.0,4375.0],[3975.0,4150.0],[3970-offset,3970+offset],[3889-offset,3889+offset],[3835-offset,3835+offset]])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    offset = 50
    lowerLine = lineList - offset

    ## cuts
    for i in range(len(lineList)):
        
        upperLine = lineWindows[i][1]
        lowerLine = lineWindows[i][0]
        cutList = np.where((wl >= lowerLine) & (wl <= upperLine))
        #cutList = np.where((wl >= lineWindows[i][0]) & (wl <= lineWindows[i][1]))
        flux = np.delete(flux,cutList)
        wl = np.delete(wl,cutList)

    ## interpolation
    #interp = sp.splrep(wl,flux,s=5,k=3)
    interp = np.interp(Owl,wl,flux)
    fluxnew = savitzky_golay(interp,103,5)
    #fluxnew = sp.splev(Owl,interp,der=0)
    #fluxnew = savitzky_golay(flux,805,5)
    #fluxsg = savitzky_golay(flux,103,5)
    #interp = sp.splrep(wl,fluxsg,s=0)
    #fluxnew = sp.interp1d(Owl,fluxnew)
    #fluxnew = sp.splev(Owl,interp,der=0)

    #fluxnew = np.interp(Owl,wl,flux)
    
    #interp = sp.InterpolatedUnivariateSpline(wl,fluxsg,k=1)
    #fluxnew = interp(Owl)
    
    #plt.plot(Owl,Oflux,alpha=0.3)
    
    #plt.plot(wl,flux,alpha=0.5)
    #plt.plot(Owl,fluxnew)
    #plt.ylim(min(Oflux),max(Oflux))
    #plt.plot(Owl,fluxnew,alpha=0.75)
    #plt.plot(Owl,testNorm)
    #plt.show()
    #plt.clf()
    
    #plt.plot(wl,flux,alpha=0.2)
    #plt.plot(Owl,fluxnew,)
    #plt.plot(Owl,Oflux,alpha=0.3)
    #plt.ylim(0,30)
    #plt.show()
    
    normalization = Oflux/fluxnew

    ####NEW STUFF
    #normalization = ((normalization - 1)*2)+1
    
    return (Owl,normalization)


def CSVGetAllVelocities(path):
    import numpy as np
    from astropy.io import fits
    import os
    import matplotlib.pyplot as plt

    c = 299792.458 #km/s
    
    Owl,Normflux = CSVNormNoPlot(path)

    #plt.plot(Owl,Normflux)
    #plt.show()

    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10

    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    

    #lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"]
    )
    for lineIndex in range(len(lineList)):
        
        #offset = 50
        #offset = 25
        #offset = 20
        #offset = 30
        
        upperLine = lineWindows[lineIndex][1]
        lowerLine = lineWindows[lineIndex][0]

        #upperLine = lineList[lineIndex] + offset
        #lowerLine = lineList[lineIndex] - offset
        
        #plt.axvline(upperLine,color='black')
        #plt.axvline(lowerLine,color="black")
        #plt.plot(Owl,Normflux)
        #plt.show()
    
        wherr = np.where((Owl >= lowerLine) & (Owl <= upperLine))
        flux = Normflux[wherr]
        
        wl = Owl[wherr]
        #wl = np.linspace(lowerLine,upperLine,len(flux))
        #plt.plot(wl,flux)
        #plt.show()
    
        vel = []
        for w in range(len(wl)):
            v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
            vel.append(v)
        if (lineIndex == 0):
            velBeta = vel
            fluxBeta = flux
            #ferrBeta = ferr
        if (lineIndex == 1):
            velGamma = vel
            fluxGamma = flux
            #ferrGamma = ferr
        if (lineIndex == 2):
            velDelta = vel
            fluxDelta = flux
            #ferrDelta = ferr

    
    vels = np.array([np.array(velBeta),np.array(velGamma),np.array(velDelta)])
    fluxes = np.array([fluxBeta,fluxGamma,fluxDelta])


    return (vels,fluxes)


def planck(wl,T):
    import numpy as np
    h = 6.626 * 10**(-27) #erg s
    c = 2.99792458 * 10**(10) #cm/s
    kb = 1.380658 * 10**(-16) #erg/K

    calc = ((2*h*c**2)/(wl**5)) * (1/(np.exp((h*c)/(wl*kb*T)) -1))
    return calc

def SDSSNormNoPlot(path):
    from astropy.io import fits
    import numpy as np
    import scipy.interpolate as sp
    import matplotlib.pyplot as plt
    from astropy.io import fits
    
    sdssData = fits.getdata(path)
    header = fits.getheader(path)
    coeff0 = header['COEFF0']
    coeff1 = header['COEFF1']
    #print sdssData['FLUX']
    flux = sdssData['FLUX']

    tmpArr = []
    for i in range(len(flux)):
        tmpArr.append(i)
    tmpArr = np.array(tmpArr)

    wl = 10**(coeff0+coeff1*tmpArr)+header['HELIO_RV']+header['SHIFT']
    Owl = 10**(coeff0+coeff1*tmpArr)+header['HELIO_RV']+header['SHIFT']
    #wl = VAC / (1.0 + 2.735182E-4 + (131.4182 / VAC**2) + (2.76249E8 / VAC**4))+header['HELIO_RV']
    #Owl = VAC / (1.0 + 2.735182E-4 + (131.4182 / VAC**2) + (2.76249E8 / VAC**4))+header['HELIO_RV']

    #wl = np.linspace(3800,9200,len(flux))
    #Owl = np.linspace(3800,9200,len(flux))
    #plt.plot(wl,flux)
    #plt.show()
    err = 1/np.array(sdssData['IVAR'])
    #print err
    Oflux = np.array(sdssData['FLUX'])
    
    #flux = flux /max(flux)
    #Oflux = Oflux / max(Oflux)
    
    #lineList =  np.array([6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397])
    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    #lineWindows = np.array([[6563.0-offset,6563.0+offset],[4800.0,4950.0],[4300.0,4375.0],[3975.0,4150.0],[3970-offset,3970+offset],[3889-offset,3889+offset],[3835-offset,3835+offset]])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    offset = 50
    lowerLine = lineList - offset
    #plt.plot(Owl,Oflux)
    ## cuts
    for i in range(len(lineList)):
        
        upperLine = lineWindows[i][1]
        lowerLine = lineWindows[i][0]
        cutList = np.where((wl >= lowerLine) & (wl <= upperLine))
        #cutList = np.where((wl >= lineWindows[i][0]) & (wl <= lineWindows[i][1]))
        flux = np.delete(flux,cutList)
        wl = np.delete(wl,cutList)
        #plt.axvline(lineList[i])
    #plt.show()
    #plt.clf()
    ## interpolation
    #interp = sp.splrep(wl,flux,s=5,k=3)
    
    interp = np.interp(Owl,wl,flux)
    fluxnew = savitzky_golay(interp,103,5)
    #fluxnew = sp.splev(Owl,interp,der=0)
    #fluxnew = savitzky_golay(flux,805,5)
    #fluxsg = savitzky_golay(flux,103,5)
    #interp = sp.splrep(wl,fluxsg,s=0)
    #fluxnew = sp.interp1d(Owl,fluxnew)
    #fluxnew = sp.splev(Owl,interp,der=0)

    #fluxnew = np.interp(Owl,wl,flux)
    
    #interp = sp.InterpolatedUnivariateSpline(wl,fluxsg,k=1)
    #fluxnew = interp(Owl)
    
    #plt.plot(Owl,Oflux,alpha=0.3)
    #plt.plot(wl,flux,alpha=0.5)
    #plt.plot(Owl,fluxnew)
    #plt.ylim(min(Oflux),max(Oflux))
    #plt.plot(Owl,fluxnew,alpha=0.75)
    #plt.plot(Owl,testNorm)
    #plt.show()
    #plt.clf()
    
    #plt.plot(wl,flux,alpha=0.2)
    #plt.plot(Owl,fluxnew,)
    #plt.plot(Owl,Oflux,alpha=0.3)
    #plt.ylim(0,30)
    #plt.show()
    
    normalization = Oflux/fluxnew
    normErrs = err / fluxnew
    ####NEW STUFF
    #normalization = ((normalization - 1)*2)+1
    #plt.plot(Owl,normalization)
    #plt.show()
    return (Owl,normalization,normErrs)


def SDSSGetAllVelocities(path):
    import numpy as np
    from astropy.io import fits
    import os
    import matplotlib.pyplot as plt

    c = 299792.458 #km/s
    
    Owl,Normflux,NormErr = SDSSNormNoPlot(path)

    #plt.plot(Owl,Normflux)
    #plt.show()

    ## Halpha, Hbeta, Hgamma, Hdelta, Hepsilon, H9, H10

    lineList = np.array([4101.734,4340.472,4861.35])
    lineWindows = np.array([[4040.0, 4190.0], [4250.0,4475.0], [4700.0,4975.0]])
    

    #lineList = np.array([4101.734,4340.472,4861.35])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4375.0], [4800.0,4950.0]])
    #lineWindows = np.array([[4060.0, 4150.0], [4300.0,4400.0], [4800.0,4950.0]])
    lineNames = np.array(["Halpha","Hbeta","Hgamma","Hdelta","Hepsilon","H9","H10"]
    )
    for lineIndex in range(len(lineList)):
        
        #offset = 50
        #offset = 25
        #offset = 20
        #offset = 30
        
        upperLine = lineWindows[lineIndex][1]
        lowerLine = lineWindows[lineIndex][0]

        #upperLine = lineList[lineIndex] + offset
        #lowerLine = lineList[lineIndex] - offset
        
        #plt.axvline(upperLine,color='black')
        #plt.axvline(lowerLine,color="black")
        #plt.plot(Owl,Normflux)
        #plt.show()
    
        wherr = np.where((Owl >= lowerLine) & (Owl <= upperLine))
        flux = Normflux[wherr]
        ferr = NormErr[wherr]
        wl = Owl[wherr]
        #wl = np.linspace(lowerLine,upperLine,len(flux))
        #plt.plot(wl,flux)
        #plt.show()
    
        vel = []
        for w in range(len(wl)):
            v = c*(lineList[lineIndex] - wl[w])/lineList[lineIndex]
            vel.append(v)
        if (lineIndex == 0):
            velBeta = vel
            fluxBeta = flux
            ferrBeta = ferr
        if (lineIndex == 1):
            velGamma = vel
            fluxGamma = flux
            ferrGamma = ferr
        if (lineIndex == 2):
            velDelta = vel
            fluxDelta = flux
            ferrDelta = ferr

    
    vels = np.array([np.array(velBeta),np.array(velGamma),np.array(velDelta)])
    fluxes = np.array([fluxBeta,fluxGamma,fluxDelta])
    ferrs = np.array([ferrBeta,ferrGamma,ferrDelta])
    
    return (vels,fluxes,ferrs)

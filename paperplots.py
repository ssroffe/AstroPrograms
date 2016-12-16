def sine(t,A,P,Phi,Gamma):
    import numpy as np
    return A*np.sin((2*np.pi)*(t/P) + Phi) + Gamma

def splitAxis(num):
    import matplotlib.pyplot as plt
    axes = plt.subplots(1,num,sharey=True,facecolor='w',figsize=(25,10))
    return  axes

def setFig():
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.rcParams.update({'figure.autolayout':'True'})
    #plt.rcParams.update({'font.size': 16})
    plt.rcParams.update({'font.size': 36})
    plt.rcParams.update({'mathtext.default':'regular'})
    plt.rcParams.update({'mathtext.fontset':'stixsans'})
    plt.rcParams.update({'axes.linewidth': 3.0})
    plt.rcParams.update({'xtick.major.size': 5})
    plt.rcParams.update({'xtick.major.width': 1.25 })
    plt.rcParams.update({'xtick.minor.size': 2.5})
    plt.rcParams.update({'xtick.minor.width': 1.25 })
    plt.rcParams.update({'ytick.major.size': 5})
    plt.rcParams.update({'ytick.major.width': 1.25 })
    plt.rcParams.update({'ytick.minor.size': 2.5})
    plt.rcParams.update({'ytick.minor.width': 1.25 })
    #plt.figure(1,figsize = [11.0, 8.5])

def TimePlot():
    import tools as tls
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    setFig()
    tls.mkdir("../../PaperPlots")

    
    lines = [line.rstrip('\n') for line in open('filelist')]

    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]
    
    tls.mkdir("../../PaperPlots/"+wdName)
    
    rvdata = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')
    
    timeArr = rvdata[:,0]
    rvArr = rvdata[:,1]
    stdArr = rvdata[:,2]

    SNArr = Signal2Noise()
    
    SNCut = 3.0
    #SNCut = 2.0
    wherrSN = np.where(SNArr >= SNCut)
    
    timeArr = timeArr[wherrSN]
    rvArr = rvArr[wherrSN]
    stdArr = stdArr[wherrSN]
    SNArr = SNArr[wherrSN]
    
    sineData = np.genfromtxt("BICFits/"+wdName+"_sineParams.csv")
    amp,period,phi,gam = sineData
    largeTime = np.linspace(np.min(timeArr)-0.2, np.max(timeArr)+0.2, 5000)
    
    sineVals = sine(largeTime, sineData[0], sineData[1], sineData[2], sineData[3])
    
    wherr1 = np.where((np.min(timeArr) <= timeArr) & (timeArr <= np.min(timeArr) + 0.5))
    #s2n = Signal2Noise()
    ## separate large time gaps into different axes
    i = 0
    wherrArr = []
    while i != len(timeArr):
        wherr = np.where((timeArr[i] <= timeArr) & (timeArr <= timeArr[i] + 0.5))
        wherrArr.append(timeArr[wherr])
        i = np.max(wherr) + 1
        
    fig, axes = splitAxis(len(wherrArr))
    #plt.xlabel("MJD [days]")
    
    ## Plotting
    off = 55950
    for i in range(len(axes)):
        axes[i].errorbar(timeArr-off,rvArr,yerr=stdArr,ls="None",marker='o',markersize=10)
        axes[i].plot(largeTime-off,sineVals,color='k')
        #ax2 = axes[i].twinx()
        #ax2.plot(timeArr,SNArr,ls="None",marker='o',color='r',markersize=10)
        #ax2.set_ylabel("Signal to Noise Ratio", color='r')
        #for t in ax2.get_yticklabels():
        #    t.set_color("red")
        if len(wherrArr[i]) == 1:
            axes[i].set_xlim(min(wherrArr[i])-off-period,max(wherrArr[i])-off+period)
            axes[i].xaxis.set_ticks([wherrArr[i]-off])
        else:
            axes[i].set_xlim(min(wherrArr[i])-off-0.005,max(wherrArr[i])-off+0.005)
            axes[i].xaxis.set_ticks(np.arange(min(wherrArr[i])-off,max(wherrArr[i])-off,0.02))
        axes[i].yaxis.tick_left()
        axes[i].ticklabel_format(useOffset=False)
        axes[i].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        #axes[i].set_ylim(-400,400)
        axes[i].set_ylim(-350,350)
        #axes[i].xaxis.get_offset_text().set_visible(False) #remove scientific notation
        axes[i].set_xlabel("MJD [days - "+str(off)+"]")
        plt.setp(axes[i].get_xticklabels(), rotation=25, horizontalalignment='right')
        #if i == int(len(axes)/2):
        #
        #    axes[i].set_title(wdName)
        if i == 0:
            axes[i].set_ylabel("RV [km/s]")

    plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_time.pdf")
    plt.savefig("../../PaperPlots/"+wdName+"/"+wdName+"_time.pdf")
    
    
    #plt.show()

def PhasePlot():
    import tools as tls
    import numpy as np
    import matplotlib.pyplot as plt
    from plot_format import plot_format
    from matplotlib.ticker import MaxNLocator
    #setFig()
    tls.mkdir("../../PaperPlots")

    lines = [line.rstrip('\n') for line in open('filelist')]

    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]

    tls.mkdir("../../PaperPlots/"+wdName)
    
    rvdata = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')
    
    timeArr = rvdata[:,0]
    rvArr = rvdata[:,1]
    stdArr = rvdata[:,2]

    SNArr = Signal2Noise()

    SNCut = 3.0
    #SNCut = 2.0
    wherrSN = np.where(SNArr >= SNCut)
    
    timeArr = timeArr[wherrSN]
    rvArr = rvArr[wherrSN]
    stdArr = stdArr[wherrSN]
    
    sineData = np.genfromtxt("BICFits/"+wdName+"_sineParams.csv")
    
    largeTime = np.linspace(np.min(timeArr), np.max(timeArr), 5000)
    
    sineVals = sine(largeTime, sineData[0], sineData[1], sineData[2], sineData[3])
    AFit,PFit,PhFit,GamFit = sineData
    
    phiDiagArr = []

    for pt in timeArr:
        PhiOff = ((2*np.pi)*(pt/PFit) + PhFit) 
        phiDiagArr.append(PhiOff)
        
        for i in range(len(phiDiagArr)):
            while (phiDiagArr[i] < 0) or (phiDiagArr[i] > (2*np.pi)):
                if (phiDiagArr[i] < 0):
                    phiDiagArr[i] = phiDiagArr[i] + (2*np.pi)
                else:
                    phiDiagArr[i] = phiDiagArr[i] - (2*np.pi)

    phiDiagArr = np.array(phiDiagArr)
    angles = np.linspace(0,2*np.pi,5000)

    yvalues = sine(phiDiagArr,AFit,2*np.pi,0,GamFit)

    residuals = rvArr - yvalues

    plot_format()
    setFig()
    plt.figure(1).add_axes((.1,.3,.8,.6))
    plt.plot(angles,sine(angles,AFit,2*np.pi,0,GamFit),'k--')
    plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
    plt.gca().xaxis.set_ticks([])
    plt.title(wdName+" Phase")
    plt.ylabel("RV [km/s]")
    plt.xlim(0,2*np.pi)
    plt.errorbar(phiDiagArr,rvArr,yerr=stdArr,linestyle='None',marker='o',markersize=10)

    ## Residuals
    plt.figure(1).add_axes((.1,.1,.8,.2))
    plt.gca().xaxis.set_ticks(np.arange(min(angles),max(angles),np.pi/2))
    plt.gca().yaxis.set_ticks([-150,0,150])
    plt.errorbar(phiDiagArr,residuals,yerr=stdArr,linestyle="None",marker="o",markersize=10)
    #plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
    plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
    plt.xlabel("Phase [rad]")
    plt.xlim(0,2*np.pi)
    #plt.ylim(-200,200)
    plt.axhline(0,linestyle="--",color="black",alpha=0.5)
    plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
    plt.savefig("../../PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
    #plt.show()

def Signal2Noise():
    from numpy import array, where, median, abs
    import tools as tls
    lines = [line.rstrip('\n') for line in open('filelist')]
    s2n = []
    for j in range(len(lines)):
        wl, flux = tls.RawSpectrum(lines[j])
        
        flux = array(flux)
        
        # Values that are exactly zero (padded) are skipped
        flux = array(flux[where(flux != 0.0)])
        n    = len(flux)
        # For spectra shorter than this, no value can be returned
        if (n>4):
            signal = median(flux)               
            noise  = 0.6052697 * median(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
            s2n.append(float(signal/noise))
    return array(s2n)

def LatexTable():
    import numpy as np
    import tools as tls
    from astropy.io import ascii

    #ascii.write(data, format='latex')
    
    lines = [line.rstrip('\n') for line in open('filelist')]
    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]

    rvdata = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')

    Objects = {"wd1235" : "J123549.89+154319.3", "wd1203" : "J120315.22+650524.4",
               "wd1140" : "J114024.02+661842.2", "wd1121" : "J112105.23+644336.4",
               "wd0907" : "J090751.78+071844.6", "wd0343" : "J034319.09+101238.0"}

    fullName = Objects[wdName]

    plusminus = "\\pm"
    timeArr = rvdata[:,0]
    rvArr = rvdata[:,1]
    stdArr = rvdata[:,2]

    nameCol = []
    rvCol = []
    for i in range(len(rvArr)):
        rvVal = "{0:.2f}".format(rvArr[i])
        stdVal = "{0:.2f}".format(stdArr[i])
        rvCol.append(str(rvVal) + " " + plusminus + " \ " + str(stdVal))
        if i == 0:
            nameCol.append(fullName)
        else:
            nameCol.append("...")
            
    ascii.write([nameCol,timeArr,rvCol],format="latex")
    
    

def PlotAll():
    TimePlot()
    PhasePlot()
    
if __name__ == '__main__':
    #TimePlot()
    #PhasePlot()
    PlotAll()
    #Signal2Noise()
    #LatexTable()

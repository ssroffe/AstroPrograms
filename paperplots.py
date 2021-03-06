### STUFF FOR PAPER

"""Get AIC Table"""
def AICTable():
    import numpy as np
    from astropy.io import ascii
    import tools as tls

    lines = [line.rstrip('\n') for line in open('filelist')]

    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]

    deltaAIC = np.genfromtxt("AICFits/"+wdName+"_deltaAIC.csv",delimiter=',')
    noOrbAIC, sineAIC, dAIC = deltaAIC

    numSpec = len(lines)

    AllObjects = {"wd0117":"J011721.34+311650.9","wd0204":"J020439.19+220724.7","wd0309":"J030941.46+005259.4",
                  "wd0332":"J033205.53+011206.7","wd0343":"J034319.09+101238.0","wd0749":"J074904.33+422420.0",
                  "wd0859":"J085921.90+043812.3","wd0907":"J090751.78+071844.6","wd0912":"J091215.43+011958.8",
                  "wd1115":"J111501.16-124217.9","wd1121":"J112105.23+644336.4","wd1137":"J113709.84+003542.9",
                  "wd1140":"J114024.02+661842.2","wd1203":"J120315.22+650524.4","wd1235":"J123549.89+154319.3",
                  "wd1306":"J130646.51+152224.9","wd1331":"J133137.06+010632.3","wd1403":"J140327.76+002119.5",
                  "wd1511":"J151132.20+451732.6","wd1521":"J152125.02+391536.5","wd1541":"J154126.42+371647.6",
                  "wd1659":"J165923.87+643809.3","wd2051":"J205118.90+031209.4","wd2229":"J222903.69+122928.6",
                  "wd2349":"J234902.80+355301.0" }

    fullname = AllObjects[wdName]

    print(fullname + " & " + str(numSpec) + " & " + str(dAIC))

"""Get Teff and Logg with errors from the Kleinman Catalog"""
def TeffLogg():
    import numpy as np
    from astropy.io import ascii
    
    Objects = {"wd2349" : "234902.80+355301.0","wd1235" : "123549.89+154319.3", "wd1203" : "120315.22+650524.4",
               "wd1140" : "114024.02+661842.3", "wd1121" : "112105.25+644336.2",
               "wd0907" : "090751.78+071844.6", "wd0343" : "034319.09+101238.0"}

    tableHead = { "SDSS":0, "Teff":28,"TeffErr":29,"logg":30,"loggErr":31 }
    #kdata = np.genfromtxt("/home/seth/research/Paperwds/KleinmanCatalog.dat")

    catalog = "/home/seth/research/Paperwds/KleinmanCatalog.dat"

    for key in Objects:
        with open(catalog) as f:
            foundObj = False
            for i, line in enumerate(f):
                splitline = line.split()
                if splitline[0] == Objects[key]:
                    objData = splitline
                    foundObj = True
                    #print splitline
                    break
            if not foundObj:
                print(key)
                raise Exception("The object is not found in the catalog")
            
        #objData = kdata[np.where(kdata[0] == Objects[key])]
        if key == 'wd2349':
            objTeff = float(objData[tableHead["Teff"]-1])
            objTeffErr = float(objData[tableHead["TeffErr"]-1])
            objLogg = float(objData[tableHead["logg"]-1])
            objLoggErr = float(objData[tableHead["loggErr"]-1])
        else:
            objTeff = float(objData[tableHead["Teff"]])
            objTeffErr = float(objData[tableHead["TeffErr"]])
            objLogg = float(objData[tableHead["logg"]])
            objLoggErr = float(objData[tableHead["loggErr"]])

        dataArr = np.array([objTeff,objTeffErr,objLogg,objLoggErr])
        print(key, dataArr)
        np.savetxt("/home/seth/research/Paperwds/"+key+"/AICFits/"+key+"_TeffLogg.csv",dataArr,delimiter=',')
    
    
"""Binary Mass function Calculation"""
def BinMassFunc():
    import tools as tls
    import numpy as np
    
    lines = [line.rstrip('\n') for line in open('filelist')]
    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]
    
    rvdata = np.genfromtxt("AICFits/"+wdName+"_rvdata.csv",delimiter=',')
    sineData = np.genfromtxt("AICFits/"+wdName+"_sineParams.csv")

    amp,period,phi,gam = sineData

    period = period * 24 * 3600 #Convert days to seconds
    amp = amp * 1000 #Convert km/s to m/s
    
    G = 6.67408 * 10**(-11) #Nm^2/kg^2
    
    f = (period*amp**3)/(2*np.pi*G) #kg

    Msf = f/(1.989 * 10**30)
    
    print("Binary mass function output: " + str(Msf) + " Msun")

    tmpArr = np.array([Msf])
    np.savetxt("AICFits/"+wdName+"_BinMassFuncVal.csv",tmpArr,delimiter=',')
    
    return f

def getCoolingModelData(path):
    import numpy as np
    from os.path import basename

    filelines = []
    Hfilelines = []
    Hefilelines = []
    base = basename(path)
    HHeFlag = False
    Cflag, COFlag, HFlag, HeFlag = False,False,False,False

    if "C_" in path:
        datastart = 5
        Cflag = True
        massVal = float(base[2] + "." + base[3])
    elif "CO_" in path:
        datastart = 0
        COFlag = True
        massVal = float(base[3] + "." + base[4])
    elif "Table_" in path:
        HHeFlag = True
        massVal = float(base[-3:])

        
    if HHeFlag:
        HData = []
        HeData = []
        with open(path) as f:
            for i, line in enumerate(f):
                if "hydrogen" in line:
                    HFlag = True
                    HIndex = i
                elif "helium" in line:
                    HFlag = False
                    HeFlag = True
                    HeIndex = i
                if HFlag and i >= HIndex + 2:
                    tmpLine = line.split()
                    Hfilelines.append(tmpLine)
                elif HeFlag and i >= HeIndex + 2:
                    tmpLine = line.split()
                    Hefilelines.append(tmpLine)
                    
        HData = np.array(Hfilelines)
        HeData = np.array(Hefilelines)
        return (HData,HeData,Cflag,COFlag,HHeFlag,massVal)
            
        
                
    else:
        with open(path) as f:
            for i, line in enumerate(f):
                if i >= datastart and len(line.split()) > 1:
                    tmpLine = line.split()
                    filelines.append(tmpLine)

        data = []
        for i in range(len(filelines)):
            try:
                if filelines[i] == filelines[i+1]:
                    del filelines[i+1]
            except:
                tmp = 1
            
        for i in range(len(filelines)):
            try:
                if i % 3 == 0:
                    newLine = filelines[i] + filelines[i+1] + filelines[i+2]
                    data.append(newLine)
            except IndexError:
                print(path)
                print(filelines[i])
                raise

    
        return (np.array(data),Cflag,COFlag,HHeFlag,massVal)
    
def CoolingModelMass():
    import tools as tls
    import numpy as np
    from plot_format import plot_format
    import matplotlib.pyplot as plt

    coolingPath = '/home/seth/research/Paperwds/coolingmodels/'
    #coolingFile = "/home/seth/research/Paperwds/coolingmodels/C_0200204"
    lines = [coolingPath+line.rstrip('\n') for line in open(coolingPath+'filelist')]
    #Objects = {"wd2349" : "234902.80+355301.0","wd1235" : "123549.89+154319.3", "wd1203" : "120315.22+650524.4",
    #           "wd1140" : "114024.02+661842.3", "wd1121" : "112105.25+644336.2",
    #           "wd0343" : "034319.09+101238.0"}
    Objects = {"wd2349" : "234902.80+355301.0","wd1235" : "123549.89+154319.3", 
               "wd1140" : "114024.02+661842.3", "wd1121" : "112105.25+644336.2",
               "wd0343" : "034319.09+101238.0"}

    
    plot_format()
    setFig()
    CTeff, Clogg = [],[]
    COTeff,COlogg = [],[]
    HTeff,Hlogg = [],[]
    HeTeff,Helogg = [],[]
    Cmass = []
    COmass = []
    HHEmass = []
    AgeArr = []
    massArr = []
    HeAgeArr = []
    
    firstCFlag = True
    firstCOFlag = True
    firstHHEFlag = True
    pnineFlag = 0
    count = 0
    flag = True
    for coolingFile in lines:

        getData = getCoolingModelData(coolingFile)
        if getData[-2]:
            HData,HeData,CFlag,COFlag,HHeFlag,massVal = getData
            HTeffData = HData[:,0]
            HloggData = HData[:,1]
            HeTeffData = np.array(HeData[:,0]).astype(np.float)
            HeloggData = np.array(HeData[:,1]).astype(np.float)
            HeAgeData = np.array(HeData[:,-1]).astype(np.float)
            
        else:
            data,Cflag,COFlag,HHeFlag,massVal = getCoolingModelData(coolingFile)
            Teff = np.array(data[:,1]).astype(np.float)
            logg = np.array(data[:,2]).astype(np.float)
            Age = np.array(data[:,4]).astype(np.float)
            
        if Cflag:
            CTeff.append(Teff)
            Clogg.append(logg)
            co = 'k'
            Cmass.append(massVal)
            #if firstCFlag:
            #    plt.plot(Teff,logg,color=co,linewidth=1.5,label="C Core")
            #    firstCFlag = False
            #else:
            #    plt.plot(Teff,logg,color=co,linewidth=1.5)
        
        elif COFlag:
            COmass.append(massVal)
            COTeff.append(Teff)
            COlogg.append(logg)
            AgeArr.append(Age)
            massArr.append(massVal)
            
            co = 'b'
            count += 1
            if firstCOFlag and massVal >= 0.45 and (count %2) == 1:
                tmp = -1
                #plt.plot(Teff,logg,color=co,linewidth=1.5,label="CO Core")
                #if not massVal == 1.0:
                #    plt.annotate(str(massVal)+"$M_\odot$",xy=(np.min(Teff),np.max(logg)+0.02),fontsize=24)
                #firstCOFlag = False
            #elif count%2 == 1:
            elif massVal >= 0.45 and massVal < 1.0 and (count % 2) == 0:
                #print coolingFile
                if pnineFlag < 2:
                    if massVal == 0.9:
                        pnineFlag += 1
                    if massVal == 0.5 and count %2 == 0:
                        plt.plot(Teff,logg,color=co,linewidth=1.5,label="CO Core")
                    elif count % 2 == 0 and not massVal == 1.0 and not massVal == 0.5:
                        if flag:
                            plt.plot(Teff,logg,color=co,linewidth=1.5)
                            plt.annotate(str(massVal)+"$M_\odot$",xy=(np.min(Teff),np.max(logg)+0.01),fontsize=22)
                        flag = True
                        if massVal == 0.9:
                            flag = False

        
        elif HHeFlag:
            HHEmass.append(massVal)
            HTeff.append(HTeffData)
            HeTeff.append(HeTeffData)
            Hlogg.append(HloggData)
            Helogg.append(HeloggData)

            HeAgeArr.append(HeAgeData)

            
            
            if firstHHEFlag and massVal <= 0.5 and massVal > 0.2:
                #plt.plot(HTeffData,HloggData,color='g',linewidth=1.5,label="H Core")
                plt.plot(HeTeffData,HeloggData,color='r',linewidth=1.5,label="He Core")
                #plt.annotate(str(massVal)+"$M_\odot$",xy=(np.min(HeTeffData),np.max(HeloggData)),fontsize=24)
                firstHHEFlag = False
                plt.annotate(str(massVal)+"$M_\odot$",xy=(np.min(HeTeffData),np.max(HeloggData)),fontsize=24)
            elif massVal <= 0.5 and massVal > 0.2:
                #plt.plot(HTeffData,HloggData,color='g',linewidth=1.5)
                plt.plot(HeTeffData,HeloggData,color='r',linewidth=1.5)
                plt.annotate(str(massVal)+"$M_\odot$",xy=(np.min(HeTeffData),np.max(HeloggData)),fontsize=24)
                #plt.annotate(str(massVal)+"$M_\odot$",xy=(np.min(HeTeffData),np.max(HeloggData)),fontsize=24)                


    CTeff,Clogg = np.array(CTeff), np.array(Clogg)
    COTeff,COlogg = np.array(COTeff), np.array(COlogg)
    HTeff,Hlogg = np.array(HTeff), np.array(Hlogg)
    HeTeff,Helogg = np.array(HeTeff), np.array(Helogg)

    ## Plot the points
    for key in Objects:
        teff,teffErr,logg,loggErr = np.genfromtxt("/home/seth/research/Paperwds/"+key+"/AICFits/"+key+"_TeffLogg.csv",delimiter=',')
        #plt.errorbar(teff,logg,xerr=teffErr,yerr=loggErr,color='green',marker='o',linewidth=2.5,markersize=11,markeredgecolor='k',markeredgewidth=1.5)
        if key == "wd1203":
            plt.annotate(key, xy=(teff-7500,logg+0.02),fontsize=22)
        elif key == 'wd0343':
            plt.annotate(key, xy=(teff+500,logg+0.035),fontsize=22)
        else:
            plt.annotate(key, xy=(teff+500,logg),fontsize=22)

    ### Age Lines ###

    ### CO Core age lines

    for i in range(len(AgeArr)):
        if AgeArr[i][0] == 0.0:
            AgeArr[i] = AgeArr[i][1:]
            COTeff[i] = COTeff[i][1:]
            COlogg[i] = COlogg[i][1:]

    maxAgeArr = AgeArr[0] ##max(AgeArr,key=len)

    allTeffArr = []
    allLoggArr = []
    for i in range(0,len(maxAgeArr)):
        TeffArr = []
        loggArr = []
        count = 0
        for j in range(len(COTeff)):
            count += 1
            if len(COTeff[j]) > i and COmass[j] >= 0.5 and (count % 2) == 0:
                TeffArr.append(COTeff[j][i])
        count = 0
        for j in range(len(COlogg)):
            count += 1
            if len(COlogg[j]) > i and COmass[j] >= 0.5 and (count % 2) == 0:
                loggArr.append(COlogg[j][i])
        allTeffArr.append(np.array(TeffArr))
        allLoggArr.append(np.array(loggArr))

        #plt.plot(TeffArr,loggArr,ls='--',color='b',alpha=0.5)
        
    allTeffArr = np.array(allTeffArr)
    allLoggArr = np.array(allLoggArr)
    count = 0
    for key in Objects:
        teff,teffErr,logg,loggErr = np.genfromtxt("/home/seth/research/Paperwds/"+key+"/AICFits/"+key+"_TeffLogg.csv",delimiter=',')
        #teff, _, logg, _ = np.genfromtxt("/home/seth/research/Paperwds/"+key+"/AICFits/"+key+"_TeffLogg.csv",delimiter=',')
        flag = True
        distArr = []
        dxArr = []
        dyArr = []
        indexArr = []
        tmpAgeArr = []
        for i in range(len(allTeffArr)):
            if not len(allTeffArr[i]) == 0:
                #allTeffArr[i] = allTeffArr[i][np.where(allLoggArr < 9.0)]
                #allLoggArr[i] = allLoggArr[i][np.where(allLoggArr < 9.0)]
                #dx,dy = dist(allTeffArr[i],teff,0,0), dist(logg,allLoggArr[i],0,0)
                dx = abs(allTeffArr[i] - teff)
                dy = abs(allLoggArr[i] - logg)
                #d = dist(allTeffArr[i],allLoggArr[i],teff,logg)
                dyI = np.argmin(dy)
                #tmpAgeArr.append(AgeArr[i])
                dxArr.append(dx[dyI])
                
                #print minD
                indexArr.append(i)
                #distArr.append(minD)
                
        

        index = indexArr[np.argmin(dxArr)]
        #index = np.where( np.array(allTeffArr) == np.array(tmpx)[np.argmin(dxArr)] )
        #print index
        alp = 0.4

        plt.plot(allTeffArr[index],allLoggArr[index],color='b',ls='--',alpha=alp)
        fs = 18
        
        if key == 'wd1203':
            plt.annotate(str("{:.2E}".format(maxAgeArr[index]))+"yr",xy=(teff-8500,logg-0.08),color='blue',fontsize=fs)
        elif key=='wd1121':
            plt.annotate(str("{:.2E}".format(maxAgeArr[index]))+"yr",xy=(teff+1200,logg-0.1),color='blue',fontsize=fs)
        elif key=='wd1235' or key=='wd2349':
            tmp = 1
        else:
            plt.annotate(str("{:.2E}".format(maxAgeArr[index]))+"yr",xy=(teff-1000,logg-0.11),color='blue',fontsize=fs)

        #plt.annotate("test",xy=(np.min(allTeffArr[index]),8.6),color='blue')
        count += 1
        plt.errorbar(teff,logg,xerr=teffErr,yerr=loggErr,color='green',marker='o',linewidth=2.5,markersize=11,markeredgecolor='k',markeredgewidth=1.5)            
        
        
        
        #plt.annotate("{:.2E}".format(maxAgeArr[index]),xy=(max(TeffArr),9))
                
        #TeffArr = [x[i] for x in COTeff]
        #loggArr = [x[i] for x in COlogg]
        #plt.plot(TeffArr,loggArr,ls='--',color='k',alpha=0.5)

    ### He Core Age Lines
    ####FUCK THIS WE'RE DOING IT MANUALLY
    for key in Objects:
        if key == 'wd1235':
            wd1235teff,teffErr,wd1235logg,loggErr = np.genfromtxt("/home/seth/research/Paperwds/"+key+"/AICFits/"+key+"_TeffLogg.csv",delimiter=',')
        elif key == 'wd2349':
            wd2349teff,teffErr,wd2349logg,loggErr = np.genfromtxt("/home/seth/research/Paperwds/"+key+"/AICFits/"+key+"_TeffLogg.csv",delimiter=',')
            
    wd1235Teffs = [18000,21000,21000,22000]
    wd1235loggs = [6.652,7.246,7.578,7.806]
    wd1235HeAge = [3.224e7]

    wd2349Age = [2.026e7]
    wd2349Teffs = [24000,25000,24000,19000]
    wd2349loggs = [7.797,7.546,7.193,6.576]

    plt.plot(wd2349Teffs,wd2349loggs,color='r',ls='--',alpha=alp)
    plt.annotate(str("{:.2E}".format(wd2349Age[0]))+"yr",xy=(wd2349teff,wd2349logg+0.1),color='r',fontsize=fs)
    
    plt.plot(wd1235Teffs,wd1235loggs,color='r',ls='--',alpha=alp)
    plt.annotate(str("{:.2E}".format(wd1235HeAge[0]))+"yr",xy=(wd1235teff,wd1235logg-0.1),color='r',fontsize=fs)
        
            
    plt.xlim(0,45000)
    plt.ylim(7,8.6)
    leg = plt.legend(loc="lower left",borderaxespad=0.,fontsize=30,frameon=True,framealpha=1)
    plt.gca().xaxis.set_ticks([10000,20000,30000,40000])
    plt.ylabel("Log(g)")
    plt.xlabel("Teff [K]")
    #plt.show()
    plt.savefig("/home/seth/research/PaperPlots/TeffLoggCombinedPlot.pdf")    
    
    
"""Sine model for fitting"""
def sine(t,A,P,Phi,Gamma):
    import numpy as np
    return A*np.sin((2*np.pi)*(t/P) + Phi) + Gamma

"""FOR VELOCITY FITTING"""
def voigtModel(x, Ldepth, Lwidth, Gdepth, Gwidth):
    import numpy as np
    return 1.0-Ldepth/(1.0 + ((x)/Lwidth)**2) - Gdepth*np.exp(-(x)**2/(2*Gwidth**2))

def lnlikeModel(p,x,y,err):
    import numpy as np
    Ldepth,Lwidth,Gdepth,Gwidth = p
    return -np.sum((y-voigtModel(x,Ldepth,Lwidth,Gdepth,Gwidth))**2/(2*err))

def lnpriorModel(p):
    import numpy as np
    Ldepth,Lwidth,Gdepth,Gwidth = p
    if 0.0 < Ldepth < 1.0 and 0.0 < Lwidth < 3000.0 and 0.0 < Gdepth < 1.0 and 0.0 < Gwidth < 600.0:
        return 0.0
    return -np.inf

def lnprobModel(p, x, y,yerr):
    import numpy as np
    lp = lnpriorModel(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlikeModel(p, x, y, yerr)

def voigt(x, Ldepth, Lwidth, Gdepth, Gwidth, RVShift):
    import numpy as np
    return 1.0-Ldepth/(1.0 + ((x-RVShift)/Lwidth)**2) - Gdepth*np.exp(-(x-RVShift)**2/(2*Gwidth**2))

"""Split the large gaps in the x-axis for the time plot"""
def splitAxis(num):
    import matplotlib.pyplot as plt
    axes = plt.subplots(1,num,sharey=True,facecolor='w',figsize=(25,10))
    return  axes

"""Make the plots look pretty"""
def setFig():
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.rcParams.update({'figure.autolayout':'True'})
    #plt.rcParams.update({'font.size': 16})
    plt.rcParams.update({'font.size': 38})
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

"""Plot orbital solution vs time"""
def TimePlot():
    import tools as tls
    from sys import platform
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    from plot_format import plot_format
    plot_format()
    setFig()
    tls.mkdir("../../PaperPlots")

    
    lines = [line.rstrip('\n') for line in open('filelist')]

    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]
    
    tls.mkdir("../../PaperPlots/"+wdName)
    
    rvdata = np.genfromtxt("AICFits/"+wdName+"_rvdata.csv",delimiter=',')
    
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
    
    sineData = np.genfromtxt("AICFits/"+wdName+"_sineParams.csv",delimiter=',')
    amp,period,phi,gam = sineData[:,0]
    
    largeTime = np.linspace(np.min(timeArr)-0.2, np.max(timeArr)+0.2, 5000)
    #sineVals = sine(largeTime, sineData[0], sineData[1], sineData[2], sineData[3])
    sineVals = sine(largeTime, amp, period,phi,gam)
    
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
    #plot_format()
    
    for i in range(len(axes)):
        axes[i].tick_params(labelsize=36)
        #axes[i].tick_params({'fontsize': 36})
        axes[i].errorbar(timeArr-off,rvArr,yerr=stdArr,ls="None",marker='o',markersize=14)
        axes[i].plot(largeTime-off,sineVals,color='k',linewidth=3.0)
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
            axes[i].xaxis.set_ticks(np.arange(min(wherrArr[i])-off,max(wherrArr[i])-off,0.04))
            
        axes[i].yaxis.tick_left()
        axes[i].ticklabel_format(useOffset=False)
        axes[i].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        #axes[i].set_ylim(-400,400)
        #axes[i].set_ylim(-350,350)
        #axes[i].xaxis.get_offset_text().set_visible(False) #remove scientific notation
        axes[i].set_xlabel("MJD [days - "+str(off)+"]")
        #if i == int(len(axes)/2):
        plt.setp(axes[i].get_xticklabels(), rotation=25, horizontalalignment='right')
        #
        #    axes[i].set_title(wdName)
        if i == 0:
            axes[i].set_ylabel("RV [km/s]")
    if platform == 'cygwin':
        #plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_time.pdf")
        plt.savefig("/cygdrive/c/Users/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_time.pdf")
    else:
        plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_time.pdf")
    plt.savefig("../../PaperPlots/"+wdName+"/"+wdName+"_time.pdf")
    
    
    #plt.show()

"""Plot phase-folded orbit"""
def PhasePlot():
    import tools as tls
    from sys import platform
    import numpy as np
    import matplotlib.pyplot as plt
    from plot_format import plot_format
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter, MultipleLocator
    #setFig()
    tls.mkdir("../../PaperPlots")

    lines = [line.rstrip('\n') for line in open('filelist')]

    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]

    tls.mkdir("../../PaperPlots/"+wdName)
    
    rvdata = np.genfromtxt("AICFits/"+wdName+"_rvdata.csv",delimiter=',')
    
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
    
    sineData = np.genfromtxt("AICFits/"+wdName+"_sineParams.csv",delimiter=',')
    sineData = sineData[:,0]
    
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

    plt.figure(1).set_size_inches(30,20)
    plt.figure(1).add_axes((.1,.3,.8,.6))
    plt.plot(angles,sine(angles,AFit,2*np.pi,0,GamFit),color='k',ls='--',linewidth=8.0)
    #plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=7,prune='lower'))
    #plt.gca().yaxis.set_ticks(np.arange(min(yvalues),max(yvalues)))
    plt.gca().yaxis.set_ticks([-200,-100,0,100,200])
    plt.gca().xaxis.set_ticks([])
    plt.title(wdName+" Phase")
    plt.ylabel("RV [km/s]")
    plt.xlim(0,2*np.pi)
    plt.errorbar(phiDiagArr,rvArr,yerr=stdArr,linestyle='None',marker='o',markersize=28)

    ## Residuals
    #plt.figure(1).add_axes((.1,.15,.8,.2))
    plt.figure(1).add_axes((.1,.15,.8,.2))
    plt.gca().xaxis.set_ticks(np.arange(min(angles),max(angles)+np.pi/2,np.pi/2))
    #print np.arange(min(angles),max(angles)+np.pi/2,np.pi/2)
    plt.gca().xaxis.set_ticklabels(['0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'])
    #plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%g $\pi$'))
    #plt.gca().xaxis.set_major_locator(MultipleLocator(base=0.5))
    plt.gca().yaxis.set_ticks([-200,0,200])
    plt.errorbar(phiDiagArr,residuals,yerr=stdArr,linestyle="None",marker="o",markersize=28)
    #plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
    #plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
    plt.xlabel("Phase [rad]")
    plt.xlim(0,2*np.pi)
    #plt.ylim(-200,200)
    plt.axhline(0,linestyle="--",color="black",alpha=0.5,linewidth=8.0)
    
    if platform == 'cygwin':
        #plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
        plt.savefig("/cygdrive/c/Users/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
    else:
        plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
    #plt.savefig("/home/seth/research/tmp/"+wdName+"_phase.png")
    plt.savefig("../../PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
    #plt.show()


"""Calculate Signal to Noise for spectra"""
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


"""make the table of rv data for Latex"""
def LatexTable():
    import numpy as np
    import tools as tls
    from astropy.io import ascii

    #ascii.write(data, format='latex')
    
    lines = [line.rstrip('\n') for line in open('filelist')]
    basename = tls.GetFileName(lines[0])
    wdName = basename[0:6]

    rvdata = np.genfromtxt("AICFits/"+wdName+"_rvdata.csv",delimiter=',')
    sineData = np.genfromtxt("AICFits/"+wdName+"_sineParams.csv",delimiter=',')
    sine = sineData[:,0]
    sineErr = sineData[:,1]



    Objects = {"wd2349" : "234902.80+355301.0","wd1235" : "J123549.89+154319.3", "wd1203" : "J120315.22+650524.4",
               "wd1140" : "J114024.02+661842.3", "wd1121" : "J112105.25+644336.2",
               "wd0907" : "J090751.78+071844.6", "wd0343" : "J034319.09+101238.0"}

    fullName = Objects[wdName]

    plusminus = "\\pm"
    timeArr = rvdata[:,0]
    rvArr = rvdata[:,1]
    stdArr = rvdata[:,2]

    ### Sine params
    sineCol = []
    nameCol = []
    row = [fullName,"$" + str("{0:.2f}".format(sine[0])) + " " + plusminus + str("{0:.2f}".format(sineErr[0]))+"$","$" + str("{0:.2E}".format(sine[1])) + " " + plusminus + str("{0:.2E}".format(sineErr[1]))+"$",
           "$" + str("{0:.2f}".format(sine[2])) + " " + plusminus + str("{0:.2f}".format(sineErr[2])) + "$","$" + str("{0:.2f}".format(sine[3])) + " " + plusminus + str("{0:.2f}".format(sineErr[3]))+ "$"]
    print("&".join(row))
    #ascii.write(row,format='latex')
        
        
    ### RVS
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
            
    #ascii.write([nameCol,timeArr,rvCol],format="latex")

"""Velocity plots"""

"""get the model Parameters to plot"""
def GetModelVelocity(Hline="gamma"):
    import tools as tls
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import os
    
    Objects = {"wd1235" : "J123549.89+154319.3", "wd1203" : "J120315.22+650524.4",
               "wd1140" : "J114024.02+661842.3", "wd1121" : "J112105.25+644336.2",
               "wd0907" : "J090751.78+071844.6", "wd0343" : "J034319.09+101238.0"}

    HLines = { "beta" : 0, "gamma" : 1, "delta" : 2 }
    Hl = HLines[Hline]
    
    lines = [line.rstrip('\n') for line in open('filelist')]
    path = lines[0]
    basename = os.path.basename(path)[:-5]
    wdName = basename[0:6]
    modelFile = "../../KoesterModels/da"+open('modelVal').read().splitlines()[0]+".dk"
    modelWl,modelFlux = tls.ModelNormNoPlot(modelFile)
    tmpWl,tmpFlux,_ = tls.NormNoPlot(lines[0])
    modelVels,modelFluxes = tls.ModelGetAllVelocities(modelFile)
    modelErrs = []
    for i in range(len(modelFluxes)):
        modelErrs.append(0.01*modelFluxes[i])

    mdim,mwalkers = 4,200
    vel = modelVels[Hl]
    flux = modelFluxes[Hl]
    err = modelErrs[Hl]
    modelSampler = tls.MCMCfit(lnprobModel,args=(np.array(vel),np.array(flux),np.array(err)),nwalkers=mwalkers,ndim=mdim,burnInSteps=16000,steps=16000)
    modelSamples = modelSampler.flatchain.reshape((-1,mdim)).T

    ld = modelSamples[0].mean()
    lw = modelSamples[1].mean()
    gd = modelSamples[2].mean()
    gw = modelSamples[3].mean()

    modelParams = np.array([ld,lw,gd,gw])
    np.savetxt("AICFits/"+wdName+"_"+Hline+"_modelParams.csv",modelParams,delimiter=',')

"""Actually plot the velocity curves"""
def PlotVelocities(Hline="gamma"):
    import os
    import tools as tls
    from plot_format import plot_format
    from sys import platform
    import numpy as np
    from matplotlib.ticker import MaxNLocator
    import matplotlib.pyplot as plt
    
    HLines = { "beta" : 0, "gamma" : 1, "delta" : 2 }
    Hl = HLines[Hline]
    lines = [line.rstrip('\n') for line in open('filelist')]
    tmppath = lines[0]
    basename = os.path.basename(tmppath)[:-5]
    wdName = basename[0:6]

    
    ld,lw,gd,gw = np.genfromtxt("BICFits/"+wdName+"_"+Hline+"_modelParams.csv",delimiter=',')
    rvdata = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')
    rvdataclone = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')
    
    timeArr = rvdata[:,0]
    rvArr = rvdata[:,1]
    stdArr = rvdata[:,2]

    rvOrd = []
    for i in range(len(rvArr)):
        rvOrd.append((rvArr[i],stdArr[i],i))

    sortArr = sorted(rvOrd)#,reverse=True)
    #print sortArr
    rvArr = [rv[0] for rv in sortArr]
    stdArr = [std[1] for std in sortArr]
    orderArr = [order[2] for order in sortArr]
    off = 1.0
    #print rvArr
    textoff = 0.1
    if len(lines) >= 10:
        halfway = int(len(lines)/2)
    else:
        halfway = 999
        
    #plot_format()
    setFig()
    plot_format()
    f, (ax1,ax2) = plt.subplots(1,2)#,sharey=True)
    for j in range(len(lines)):
        path = lines[orderArr[j]]
        c = 299792.458
        vels,fluxes,ferrs = tls.GetAllVelocities(path)
        vels = vels[Hl]
        fluxes = fluxes[Hl]
        ferrs = ferrs[Hl]
        stdalph = 0.6
        modelLW = 1.0
        
        if j <= halfway:
            k = j
            ax2.step(vels,fluxes+k*off,where='mid',linewidth=1.0,color='k')
            ax2.plot(vels,voigt(vels,ld,lw,gd,gw,rvArr[j])+k*off,color='r',linewidth=modelLW,label='Best fit')
            ax2.plot(vels,voigt(vels,ld,lw,gd,gw,rvArr[j]+stdArr[j])+k*off,color='cyan',linewidth=modelLW,alpha=stdalph,label='Best fit $\pm$ std')
            ax2.plot(vels,voigt(vels,ld,lw,gd,gw,rvArr[j]-stdArr[j])+k*off,color='cyan',linewidth=modelLW,alpha=stdalph)
            ax2.axvline(0,color='k',ls='--')
            ax2.text(-1450,np.min(fluxes)+k*off-textoff,"Spec. "+str(orderArr[j]),fontsize=10)
            plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
            ax2.set_ylabel("Normalized Flux + offset")
            ax2.set_xlabel("Velocity [km s$^{-1}$]")
            ax2.set_xlim(-1500,1500)
            plt.gca().xaxis.set_ticks([-1000,0,1000])
            ax2.set_ylim(min(fluxes-0.2),max(fluxes)+k*off+0.2)
            #ax1.invert_yaxis()
        else:
            k = j - halfway
            ax1.step(vels,fluxes+k*off,where='mid',linewidth=1.0,color='k')
            ax1.plot(vels,voigt(vels,ld,lw,gd,gw,rvArr[j])+k*off,color='r',linewidth=modelLW,label='Best fit')
            ax1.plot(vels,voigt(vels,ld,lw,gd,gw,rvArr[j]+stdArr[j])+k*off,color='cyan',linewidth=modelLW,alpha=stdalph,label='Best fit $\pm$ std')
            ax1.plot(vels,voigt(vels,ld,lw,gd,gw,rvArr[j]-stdArr[j])+k*off,color='cyan',linewidth=modelLW,alpha=stdalph)
            ax1.text(-1450,np.min(fluxes)+k*off-textoff,"Spec. "+str(orderArr[j]),fontsize=10)
            ax1.axvline(0,color='k',ls='--')
            ax1.set_xlabel("Velocity [km s$^{-1}$]")
            ax1.set_ylabel("Normalized Flux + offset")
            ax1.set_xlim(-1500,1500)
            plt.gca().xaxis.set_ticks([-1000,0,1000])
            ax1.set_ylim(min(fluxes)+0.2,max(fluxes)+k*off)
    #plt.show()
    if platform == 'cygwin':
        #plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
        plt.savefig("/cygdrive/c/Users/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_velocity.pdf")
    else:
        plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")

    plt.savefig("../../PaperPlots/"+wdName+"/"+wdName+"_velocity.pdf")

def PlotOneVelocity(Hline='gamma'):
    import os
    import tools as tls
    from plot_format import plot_format
    from sys import platform
    import numpy as np
    from matplotlib.ticker import MaxNLocator, MultipleLocator
    import matplotlib.pyplot as plt
    
    HLines = { "beta" : 0, "gamma" : 1, "delta" : 2 }
    Hl = HLines[Hline]
    lines = [line.rstrip('\n') for line in open('filelist')]
    tmppath = lines[0]
    basename = os.path.basename(tmppath)[:-5]
    wdName = basename[0:6]

    minorLocator = MultipleLocator(100)
    
    ld,lw,gd,gw = np.genfromtxt("BICFits/"+wdName+"_"+Hline+"_modelParams.csv",delimiter=',')
    rvdata = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')
    rvdataclone = np.genfromtxt("BICFits/"+wdName+"_rvdata.csv",delimiter=',')
    
    timeArr = rvdata[:,0]
    rvArr = rvdata[:,1]
    stdArr = rvdata[:,2]    

    ind = np.argmax(rvArr)

    path = lines[ind]
    vels,fluxes,ferrs = tls.GetAllVelocities(path)
    vels = vels[Hl]
    fluxes = fluxes[Hl]
    ferrs = ferrs[Hl]
    rv = rvArr[ind]
    std = stdArr[ind]
    modelLW = 3.0
    vpstd = voigt(vels,ld,lw,gd,gw,rv+std)
    vnstd = voigt(vels,ld,lw,gd,gw,rv-std)
    pstd = rv+std
    nstd = rv-std
    stdalph = 0.6

    setFig()
    plot_format()
    setFig()
    plt.step(vels,fluxes,where='mid',linewidth=2.0,color='k')
    plt.plot(vels,voigt(vels,ld,lw,gd,gw,rv),color='r',linewidth=modelLW,label='Best Fit')
    plt.plot(vels,vpstd,color='cyan',linewidth=modelLW,alpha=stdalph,label='Best fit $\pm$ std')
    plt.plot(vels,vnstd,color='cyan',linewidth=modelLW,alpha=stdalph)

    plt.axvspan(nstd,pstd,color='grey',alpha=0.5)
    plt.axvline(0,color='k',ls='--',linewidth=2.0)

    plt.ylabel("Normalized Flux + offset")
    plt.xlabel("Velocity [km s$^{-1}$]")
    plt.xlim(-1500,1500)
    plt.gca().xaxis.set_ticks([-1000,-500,0,500,1000])
    plt.gca().xaxis.set_minor_locator(minorLocator)
    plt.xticks(rotation=25)

    plt.tick_params(which='both',width=1.5)
    plt.tick_params(which='minor',length=7)
    plt.ylim(min(fluxes-0.2),max(fluxes)+0.2)
    
    if platform == 'cygwin':
        #plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")
        plt.savefig("/cygdrive/c/Users/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_velocity.pdf")
    else:
        plt.savefig("/home/seth/Dropbox/astro_research/PaperPlots/"+wdName+"/"+wdName+"_phase.pdf")

    plt.savefig("../../PaperPlots/"+wdName+"/"+wdName+"_velocity.pdf")

    
    
def PlotAll():
    TimePlot()
    PhasePlot()
    PlotVelocities()

def makeVelocityCurves():
    GetModelVelocity()
    PlotVelocities()

def dist(x,y,x2,y2):
    import numpy as np
    d = np.sqrt((x-x2)**2 + (y-y2)**2)
    return d
    
if __name__ == '__main__':
    TimePlot()
    PhasePlot()
    #PlotAll()
    #LatexTable()
    #BinMassFunc()
    #GetModelVelocity()
    #PlotOneVelocity()
    #PlotVelocities()
    #CoolingModelMass()
    #TeffLogg()
    #AICTable()

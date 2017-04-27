import numpy as np
import matplotlib.pyplot as plt
import tools as tls
import os

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


lines = [line.rstrip('\n') for line in open('filelist')]
basename=os.path.basename(lines[0])[:-5]
wdName = basename[0:6]

sineData = np.genfromtxt("AICFits/"+wdName+"_sineParams.csv",delimiter=',')
sineParams = sineData[:,0]

rvData = np.genfromtxt("AICFits/"+wdName+"_rvdata.csv",delimiter=',')
timeArr = rvData[:,0]
rvArr = rvData[:,1]
stdArr = rvData[:,2]

wgtAvg = (np.sum(rvArr * stdArr**(-2))) / np.sum(stdArr**(-2))
wgtStd = 1 / (np.sum(stdArr**(-2)))
mparam = wgtAvg

noOrbParams = (mparam)
noOrbk = 1
#noOrbBIC = -2*lnlikeNoOrbit(mparam,timeArr,rvArr,stdArr)+noOrbk*np.log(len(timeArr))
noOrbAIC = -2*lnlikeNoOrbit(mparam,timeArr,rvArr,stdArr)+2*noOrbk + ( (2*noOrbk*(noOrbk+1)) / (len(timeArr) - noOrbk - 1))

#sineParams = (Afit,Pfit,Phfit,Gfit)
sinek = 4
#sineBIC = -2*lnlikeSine(sineParams,timeArr,rvArr,stdArr)+sinek*np.log(len(timeArr))
sineAIC = -2*lnlikeSine(sineParams,timeArr,rvArr,stdArr)+2*sinek + ( (2*sinek*(sinek+1)) / (len(timeArr) - sinek - 1))

deltaAIC = noOrbAIC - sineAIC

#print(noOrbAIC,sineAIC,deltaAIC)

#bicFile = open("AICFits/"+wdName+"_BICCalc.txt",'w')
#bicFile.write("Orbit eqn: v(t) = {0:.3f}*sin(2*pi*(t/{1:.3f}) + {2:.3f}) + {3:.3f}\n".format(Afit,Pfit,Phfit,Gfit))
#bicFile.write("No Orbit eqn: v(t) = {0:.3f}\n".format(float(mparam)))
#bicFile.write("No Orbit BIC = {0:.3f}\n".format(float(noOrbBIC)))
#bicFile.write("Sine AIC = {0:.3f}\n".format(float(sineBIC)))
#bicFile.write("Delta AIC = noOrbBIC - sineBIC = {0:.3f}".format(float(deltaBIC)))
#bicFile.close()
deltaAICArr = np.array([noOrbAIC, sineAIC, deltaAIC])
np.savetxt("AICFits/"+wdName+"_deltaAIC.csv",deltaAICArr,delimiter=',')

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

def plot_format():
    #Plots
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

def label(x,y,title):
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)

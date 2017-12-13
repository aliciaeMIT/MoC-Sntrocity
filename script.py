#plotting data analysis

from numpy import genfromtxt
import numpy as np
import matplotlib
import math
from plotter import mkdir_p
import time

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lines = True
nolines = False

gen_time = 'data_plots/'+ time.strftime("%Y-%m-%d_%H-%M-%S")
mkdir_p(gen_time)

def readMOCdata():
    MOCdata = genfromtxt('MOCdata.txt', delimiter='\t')

    unknowns = MOCdata[:,0]
    fm_ratio = MOCdata[:,1]
    corner_normalized = MOCdata[:,2]
    corner_ratio = MOCdata[:,3]
    fmod_error_ratio = MOCdata[:,4]
    ln_corner = MOCdata[:,5]
    ln_ratio = MOCdata[:,6]
    ln_unknowns = MOCdata[:,7]
    return unknowns, corner_ratio, fmod_error_ratio

def readSNdata():
    SNdata = genfromtxt('sndata.txt', delimiter='\t')
    unknowns = SNdata[:,0]
    fm_ratio = SNdata[:,1]
    corner_normalized = SNdata[:,2]
    corner_ratio = SNdata[:,3]
    fmod_error_ratio = SNdata[:,4]
    ln_corner = SNdata[:,5]
    ln_ratio = SNdata[:,6]
    ln_unknowns = SNdata[:,7]
    return unknowns, corner_ratio, fmod_error_ratio


def plotNoLines():
    fig1 = plt.figure(1)
    plt.loglog(newx, yvals, 'r.', label='MOC')
    plt.loglog(newx2, y2vals, 'b*', label='SN')
    #plt.ylabel('abs((new-old/new))')
    #plt.ylabel(r"$\left|\frac{\phi_{new} - \phi_{old}}{\phi_{new}}\right|$")
    plt.ylabel(r"$\mathbf{\varepsilon}$", fontsize = 24, rotation='horizontal')
    plt.xlabel(r"$\frac{1}{n}$", fontsize = 16)
    plt.title('Fuel corner flux ratio error vs 1/\# unknowns')
    plt.legend()
    plt.savefig(gen_time + '/corner_flux_error_convergence.png', dpi=1000)


    #plt.close()
    fig2 = plt.figure(2)
    plt.loglog(newx, y3vals, 'r.', label='MOC')
    plt.loglog(newx2, y4vals, 'b*', label='SN')
    #plt.ylabel(r'\textbf{abs((new-old/new))}')
    #plt.xlabel('1/ \# unknowns')
    plt.ylabel(r"$\mathbf{\varepsilon}$", fontsize = 24, rotation='horizontal')
    plt.xlabel(r"$\frac{1}{n}$", fontsize = 16)
    plt.title('Fuel-to-moderator average flux error vs 1/\# unknowns')
    plt.legend()
    plt.savefig(gen_time + '/fmod_ratio_error_convergence.png', dpi=1000)

def plotWithLines():
    fig1 = plt.figure(1)
    logx=[]
    logy=[]
    for x in newx:
        logx.append(math.log(x))
    for y in yvals:
        logy.append(math.log(y))

    fit = np.polyfit(logx, logy, 1)
    line =[]
    for lx in logx:
        line.append(np.exp(fit[0]*lx + fit[1]))
    plt.loglog(newx, yvals, 'r.', newx, line, 'r--', linewidth=0.5, label='MOC')
    print "corner ratio MOC fit:\t m \t %g\tc \t%g" %(fit[0], fit[1])

    logx2 = []
    logy2 = []
    for x in newx2:
        logx2.append(math.log(x))
    for y in y2vals:
        logy2.append(math.log(y))

    fit2 = np.polyfit(logx2, logy2, 1)

    line2 = []
    for lx in logx2:
        line2.append(np.exp(fit2[0] * lx + fit2[1]))
    plt.loglog(newx2, y2vals, 'b*', newx2, line2, 'b--', linewidth=0.5, label='SN')
    print "corner ratio SN fit:\t m \t%g\t c \t%g" %(fit2[0], fit2[1])

    plt.legend()
    plt.ylabel(r"$\mathbf{\varepsilon}$", fontsize=24, rotation='horizontal')
    plt.xlabel(r"$\frac{1}{n}$", fontsize=16)
    #plt.xlabel(r"\# unknowns", fontsize=16)
    plt.title('Fuel corner flux ratio error vs 1/\# unknowns')
    #plt.title('Fuel corner flux ratio error vs \# unknowns')
    plt.savefig(gen_time + '/corner_flux_error_convergence.png', dpi=1000)

    #*****************************************************#
    fig2 = plt.figure(2)
    #plt.loglog(newx, y3vals, 'r.', label='MOC')
    #plt.loglog(newx2, y4vals, 'b*', label='SN')

    logx1 = []
    logy1 = []
    for x in newx:
        logx1.append(math.log(x))
    for y in y3vals:
        logy1.append(math.log(y))

    fit1 = np.polyfit(logx1, logy1, 1)
    line1 = []
    for lx in logx1:
        line1.append(np.exp(fit1[0] * lx + fit1[1]))
    plt.loglog(newx, y3vals, 'r.', newx, line1, 'r--', linewidth=0.5, label='MOC')
    print "flux ratio MOC fit:\t m \t%g\t c \t%g" % (fit1[0], fit1[1])


    logx3 = []
    logy3 = []
    for x in newx2:
        logx3.append(math.log(x))
    for y in y4vals:
        logy3.append(math.log(y))

    fit3 = np.polyfit(logx3, logy3, 1)

    line3 = []
    for lx in logx3:
        line3.append(np.exp(fit3[0] * lx + fit3[1]))
    plt.loglog(newx2, y4vals, 'b*', newx2, line3, 'b--', linewidth=0.5, label='SN')
    print "flux ratio SN fit:\t m \t%g\t c \t%g" % (fit3[0], fit3[1])


    plt.ylabel(r"$\mathbf{\varepsilon}$", fontsize=24, rotation='horizontal')
    #plt.xlabel(r"\# unknowns", fontsize=16)
    plt.xlabel(r"$\frac{1}{n}$", fontsize=16)
    #plt.title('Fuel-to-moderator average flux error vs \# unknowns')
    plt.title('Fuel-to-moderator average flux error vs 1/\# unknowns')
    plt.legend()
    plt.savefig(gen_time + '/fmod_ratio_error_convergence.png', dpi=1000)

############################################################################################
############################################################################################
############################################################################################

newx = []
newx2 = []
xvals, yvals, y3vals = readMOCdata()
x2vals, y2vals, y4vals = readSNdata()
for val in xvals:
    newx.append(1 / val)
    #newx.append(val)
for val in x2vals:
    newx2.append(1 / val)
    #newx2.append(val)

if nolines:
    plotNoLines()
if lines:
    plotWithLines()
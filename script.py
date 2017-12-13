#plotting data analysis

from numpy import genfromtxt
import matplotlib.pyplot as plt

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
    return unknowns, corner_ratio

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
    return unknowns, corner_ratio

newx=[]
newx2=[]
xvals, yvals = readMOCdata()
x2vals, y2vals = readSNdata()
for val in xvals:
    newx.append(1/val)
for val in x2vals:
    newx2.append(1/val)
plt.figure()
plt.loglog(newx, yvals, 'o')
plt.loglog(newx2, y2vals, '*')
#plt.ylabel('Scalar Flux')
#plt.xlabel('X node # across centerline)
#plt.title('Horizontal centerline flux')
#plt.savefig('horiz_flux_center.png', dpi=1000)
plt.show()
#plt.close()
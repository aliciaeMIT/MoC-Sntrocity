#Alicia M. Elliott, 22.212 Fall 2017
#Method of Characteristics solver, with SN-like geometry

from initializetracks import InitializeTracks

from cellflux import MethodOfCharacteristics, ConvergenceTest
from math import pi
import geometry as geom
import plotter
import time


###################################
########## PROBLEM SETUP ##########
###################################
pitch = 1.6
fwidth = 0.8                    #fuel width/height
num_azim = 24                    #number of azimuthal angles desired
t = 0.05                         #track spacing desired, cm
spacing = 0.01                   #mesh spacing
n_p = 3                         #number of polar divisions; can be 2 or 3
num_iter_max = 100              #maximum number of iterations on flux
tol = 1e-3                      #tolerance for convergence (using L2 Engineering Norm)
fuelgeom = 'square'

h = pitch                       #height of pincell
w = pitch                       #width of pincell
r = fwidth/2                    #fuel pin effective radius (half width of square fuel pin)
#########################################
########## MATERIAL PROPERTIES ##########
#########################################

q_fuel = 10/(4 * pi)             #constant isotropic source in fuel
q_mod = 0                       #no source in moderator
update_source = False
"""
ndensity_fuel = 2.2e22          #atoms/cc (UO2)
ndensity_mod = 1.0e21           #at/cc (H2O)
sigma_a = 1e-24                 #moderator absorption cross section (cm^2)
sigma_r = (11.4 + 8)* 1e-24     #fuel absorption cross section (cm^2)
"""

#########################################
########## MATERIAL PROPERTIES ##########
#########################################
density_uo2 = 10.4                          #g/cc
density_h2o = 0.7                           #g/cc
molwt_u238 = 238                            #g/mol
molwt_o = 16                                #g/mol
molwt_h2 = 2                                #g/mol

molwt_uo2 = molwt_u238 + 2 * molwt_o
molwt_h2o = molwt_h2 + molwt_o

N_A = 6.022e23                              #1/mol; Avogadro's number

num_density_uo2 = density_uo2 * (N_A / molwt_uo2) #at/cc
num_density_h2o = density_h2o * (N_A / molwt_h2o) #at/cc

#microscopic XS, in barns
xs_scatter_238 = 11.29
xs_scatter_o = 3.888
xs_scatter_h = 20.47
xs_absorption_h = 1.0

barns = 1e-24 #cm; conversion factor

#macroscopic XS, cm^-1
sigma_fuel_scatter = (xs_scatter_238 + 2 * xs_scatter_o) * num_density_uo2 * barns
sigma_mod_scatter = (xs_scatter_o + 2 * xs_scatter_h) * num_density_h2o * barns

sigma_fuel_abs = 0.0
sigma_mod_abs = (2 * xs_absorption_h) * num_density_h2o * barns

sigma_fuel_tot = sigma_fuel_scatter + sigma_fuel_abs
sigma_mod_tot = sigma_mod_scatter + sigma_mod_abs

#########################################
############ RESULTS STORAGE ############
#########################################

#create directory to store plots, results in
timestr = time.strftime("%Y-%m-%d_%H-%M")
pathname = 'plots/' + timestr
plotter.mkdir_p(pathname)
savepath = pathname
resultsfile = pathname + '/' + timestr + '_results'

f = open('%s.txt' %resultsfile, 'w+')

f.write("********PROBLEM SETUP********\n")
f.write("cell pitch \t %g\nfuel width \t %g\nfuel source \t %g\nmod source \t %g\n\n" %(pitch, fwidth, q_fuel, q_mod))
f.write("converge tol \t %g\nnum_azim \t %g\nnum_polar \t %g\n" %(tol, num_azim, n_p))
f.write("track spacing \t %g\nmesh spacing \t %g\n" %(t, spacing))
f.write("fuel total xs \t %g\nfuel scatter \t %g\nfuel absorp \t %g\n"
        "mod total xs \t %g\nmod scatter \t %g\nmod absorp \t%g\n"
        "*****************************\n\n" %(sigma_fuel_tot, sigma_fuel_scatter, sigma_fuel_abs, sigma_mod_tot, sigma_mod_scatter,sigma_mod_abs))


###############################################
########## SETUP FLAT SOURCE REGIONS ##########
###############################################
#set material objects
fuelmat = geom.Material('fuel', q_fuel, sigma_fuel_tot, sigma_fuel_scatter, sigma_fuel_abs)
moderator = geom.Material('moderator', q_mod, sigma_mod_tot, sigma_mod_scatter, sigma_mod_abs)


fuel = geom.FlatSourceRegion(q_fuel, sigma_fuel_tot)
mod = geom.FlatSourceRegion(q_mod, sigma_mod_tot)
fsr = [fuel, mod]


# setup mesh cells
mesh = geom.Geometry(pitch, spacing, fwidth, fuelmat, moderator)
mesh.setMesh()

cell_width = mesh.getWidth(pitch, spacing)
fuel_width = mesh.getWidth(fwidth, spacing)
plot_cells = mesh.getPlotCells(cell_width, fuel_width)
plotter.plotMaterial(mesh, spacing, plot_cells, savepath)

#####################################
########## GENERATE TRACKS ##########
#####################################
check = ConvergenceTest()
setup = InitializeTracks(num_azim, t, w, h, n_p, r, fsr, fuelgeom)
setup.getTrackParams()
setup.makeTracks()
setup.findAllTrackCellIntersect(mesh.cells, spacing)
f.write("\nTotal number of segments \t %g\n\n" % (setup.tot_num_segments))
print "\nTotal number of segments \t %g\n\n" % (setup.tot_num_segments)
"""
setup.plotCellSegments(spacing, savepath)

setup.getAngularQuadrature()
setup.getPolarWeight()
#if fuelgeom == 'square':
#    setup.findIntersection()
#else:
#    setup.findCircleIntersection()
setup.reflectRays()
setup.getFSRVolumes(fuel, mod, mesh)
#setup.getTrackLinkCoords()
"""
##############################
########## PLOTTING ##########
##############################

#setup.plotTracks(savepath)
#setup.plotTrackLinking(savepath)
#setup.plotSegments()

######################################
########## SOLVE FOR FLUXES ##########
######################################

"""
flux = MethodOfCharacteristics(sigma_fuel_tot, sigma_mod_tot, fsr, setup, check, mesh)
flux.solveFlux(num_iter_max, tol, update_source)
f.write("\nConverged in %d iterations! \nAvg fuel flux\t %f \nAvg mod flux\t %f\nAverage Flux\t %f \nFlux ratio\t %f\n\n"
        % (flux.results[0], flux.results[1], flux.results[2], flux.results[3], flux.results[4]))

midpt = mesh.n_cells/2 - 1
plotter.plotScalarFlux(mesh, 0, mesh.mesh, 0, savepath)
plotter.plotCenterFlux(mesh, mesh.cells, midpt, 0, 0, savepath)
plotter.plotCenterFluxY(mesh, mesh.cells, midpt, 0, 0, savepath)
"""
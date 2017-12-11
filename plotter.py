
from PIL import Image
from PIL import ImageDraw
from math import *
import matplotlib.pyplot as plt
from errno import EEXIST
from os import makedirs,path
"""
Plotting tool for SN/MOC to plot mesh, scalar flux

plotMaterial, plotScalarFlux methods are adapted from Samuel Shaner's Sn solver, DiscOrd:
https://github.com/samuelshaner/DiscOrd/

"""

def plotCenterFlux(mesh, solved, j, iter, order, savepath):
    ivals = []
    cellfluxes = []
    fuelbounds = [(mesh.n_mod-0.5), (mesh.n_fuel + mesh.n_mod - 0.5)]

    for xc in fuelbounds:
        plt.axvline(x=xc, color='k', linestyle = '--')

    for i in range(mesh.n_cells):
        ivals.append(i)
        cellfluxes.append(solved[i][j].flux)
    plt.plot(ivals, cellfluxes)
    plt.ylabel('Scalar Flux')
    plt.xlabel('X node # across centerline (Y node num ' + str(j) + ' )')
    plt.title('Horizontal centerline flux')
    plt.savefig(savepath + '/horiz_flux_center.png')
    plt.close()

def plotCenterFluxY(mesh, solved, i, iter, order, savepath):
    jvals = []
    cellfluxes = []
    fuelbounds = [(mesh.n_mod - 0.5), (mesh.n_fuel + mesh.n_mod - 0.5)]
    #0.5 is subtracted to account for boundary being between cells for fuel and mod indices

    for xc in fuelbounds:
        plt.axvline(x=xc, color='k', linestyle = '--')

    for j in range(mesh.n_cells):
        jvals.append(j)
        cellfluxes.append(solved[i][j].flux)
    plt.plot(jvals, cellfluxes)
    plt.ylabel('Scalar Flux')
    plt.xlabel('Y node # across centerline (X node num ' + str(i) + ' )')
    plt.title('Vertical centerline flux')
    plt.savefig(savepath + '/vert_flux_center.png')
    plt.close()

def plotMaterial(mesh, spacing, plot_cells, savepath):


    bit_size = 500.0

    # create image
    bit_length = round(bit_size / mesh.n_cells)
    size = int(bit_length * mesh.n_cells)
    img = Image.new('RGB', (size,size), 'white')
    draw = ImageDraw.Draw(img)

    # draw cell interior
    for i in range(mesh.n_cells):
        for j in range(mesh.n_cells):
            cell = mesh.cells[i][j]
            # fuel red; moderator blue
            if cell.region == 'fuel':
                draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (255,0,0))
            else:
                draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (0,0,255))

            # find cells where angular flux needs to be plotted; make them white
            #if mesh.cells[i][j].id in plot_cells:
            #    draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (255,255,255))

    # draw horizontal grid lines
    for j in range(1,mesh.n_cells):
        draw.line((0, j*bit_length, size,j*bit_length), fill=(0,0,0), width=1)

    # draw vertical grid lines
    for i in range(1,mesh.n_cells):
        draw.line((i*bit_length, 0, i*bit_length, size), fill=(0,0,0), width=1)

    # save image
    img.save(savepath + '/material_' + str(spacing)[2:] + '.png')
    #img.close()

def plotScalarFlux(mesh, order, spacing, iteration, savepath):

    # create image
    bit_length = round(500.0 / mesh.n_cells)
    size = int(bit_length * mesh.n_cells)
    img = Image.new('RGB', (size,size), 'white')
    draw = ImageDraw.Draw(img)

    # find max and min flux
    max_flux = mesh.cells[0][0].flux
    min_flux = mesh.cells[0][0].flux
    for i in range(mesh.n_cells):
        for j in range(mesh.n_cells):
            max_flux = max(mesh.cells[i][j].flux, max_flux)
            min_flux = min(mesh.cells[i][j].flux, min_flux)

    flux_range = max_flux - min_flux

    # draw flux map
    for i in range(mesh.n_cells):
        for j in range(mesh.n_cells):
            cell = mesh.cells[i][j]


            # get color
            if ((cell.flux-min_flux) / flux_range <= 1.0/3.0):
                red = 0.0
                green = 3.0 * (cell.flux-min_flux) / flux_range
                blue = 1.0
            elif ((cell.flux-min_flux) / flux_range <= 2.0/3.0):
                red = 3.0 * (cell.flux-min_flux) / flux_range - 1.0
                green = 1.0
                blue = -3.0 * (cell.flux-min_flux) / flux_range + 2.0
            else:
                red = 1.0
                green = -3.0 * (cell.flux-min_flux) / flux_range + 3.0
                blue = 0.0

            # convert color to RGB triplet
            red = int(255*red)
            green = int(255*green)
            blue = int(255*blue)

            # draw pin and pin power
            draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (red,green,blue))

    # save image
    img.save(savepath + '/flux_' + str(spacing)[2:] + '_' + str(int(floor(order/10))) + str(order % 10) + '_' + str(int(floor(iteration/10))) + str(iteration % 10) + '.png')


def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

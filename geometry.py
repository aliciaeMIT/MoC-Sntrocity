from __future__ import division
import numpy as np
import math


class Geometry(object):
    def __init__(self, pitch, mesh_spacing, fuel_width, fuel, moderator):

        self.width = pitch
        self.mesh = mesh_spacing
        self.fw = fuel_width
        self.fuel = fuel
        self.moderator = moderator
        self.fuel_area = fuel_width ** 2
        self.mod_area = pitch ** 2 - fuel_width ** 2

    def setMesh(self):
        self.n_cells = int(math.ceil(self.width / self.mesh)) #number of cells in x, y
        n_cells = self.n_cells
        self.cells = np.zeros((n_cells, n_cells), dtype=object)

        self.n_fuel = int(math.ceil(self.fw / self.mesh))
        self.n_mod = int((self.n_cells - self.n_fuel) / 2)
        if not (2*self.n_mod + self.n_fuel) == self.n_cells:
            print "Check number of cells calculation"

        #setup 2d array of cells
        for j in range(n_cells):
            for i in range(n_cells):
                self.cells[i][j] = Cell()
                cell = self.cells[i][j]
                xl = self.mesh * i
                xr = xl + self.mesh
                yb = self.mesh * j
                yt = yb + self.mesh
                cell.getBoundaries(xl, xr, yt, yb)

                if (i >= self.n_mod and j >= self.n_mod) and \
                        ((i < (self.n_mod + self.n_fuel) and j < (self.n_mod + self.n_fuel))):
                    cell.region = 'fuel'
                    cell.getMaterial(self.fuel)

                    #print "set cell %d, %d to %s " % (i, j, cell.region)
                else:
                    cell.getMaterial(self.moderator)

                    #print "set cell %d, %d to %s " % (i, j, cell.region)

    def getWidth(self, width, spacing):
        return int(width / spacing)

    def getPlotCells(self, cell_width, fuel_width):
        plotcells = [int((cell_width / 2.0) * cell_width + cell_width / 2.0),
                  int((cell_width / 2.0 + fuel_width / 2 - 1) * cell_width + cell_width / 2.0) + fuel_width / 2 - 1,
                  int(cell_width * cell_width - 1),
                  int((cell_width / 2.0) * cell_width + cell_width / 2.0) + fuel_width / 2 - 1,
                  int((cell_width / 2.0 + 1) * cell_width - 1)]
        return plotcells


class Cell(object):
    def __init__(self):

        self.region = 'moderator'
        self.volume = 0
        self.area = 0
        self.flux = 0
        #self.angular = np.zeros((4, 200))
        self.source = 0
        #self.avg_angular = np.zeros(200)
        self.segments = []
        self.corr = 1

    def getMaterial(self, material):
        self.material = material

    def getBoundaries(self, xl, xr, yt, yb):
        self.xl = xl
        self.xr = xr
        self.yt = yt
        self.yb = yb

class FlatSourceRegion(object):
    def __init__(self, source, xs):
        """
        This class generates the source for each flat source
        """
        self.region = 'moderator'
        self.volume = 0
        self.area = 0
        self.flux = 0
        self.source = source
        self.sigma = xs
        # self.segments = []

    def getMaterial(self, material):
        self.material = material



class Material(object):
    def __init__(self, region, q, xs, scatter, absorption):
        self.name = region
        self.q = q
        self.xs = xs  # total xs
        self.scatter = scatter
        self.absorption = absorption
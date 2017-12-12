
import math


class MethodOfCharacteristics(object):
    """
    class for calculating angular and scalar flux
    """
    def __init__(self, sigma_t_fuel, sigma_t_mod, regions, setup, check, mesh):

        self.sigma_t_fuel = sigma_t_fuel
        self.sigma_t_mod = sigma_t_mod
        self.regions = regions
        self.setup = setup
        self.tracks = setup.tracks
        self.check = check
        self.exponential = []
        self.mesh = mesh
        self.cells = mesh.cells
        #get n_p by self.setup.n_p, numazim2 by self.setup.num_azim2
        #get segments in a track by self.tracks.segments
        #get segments in a region by self.regions.segments

    def exponentialTerm(self, s, p): #calculates exponential term for a given segment and polar angle

        ans = math.exp(-1 * s.cell.material.xs * s.length / self.setup.sintheta_p[p])
        return ans

    def preCalculate(self): #precalculates exponential terms for each segment and polar angle
        for i in range(self.setup.num_azim2):
            for track in self.tracks[i]:
                for s in track.segments:
                    for p in range(self.setup.n_p):
                        try:
                            s.exponential.append(self.exponentialTerm(s, p))
                        except(UnboundLocalError):
                            print "segment in cell %d, %d" %(s.cellid_i, s.cellid_j)
                            raise



    def angularFlux(self, flux_in, s, p): #s = segment
        q_seg = s.cell.material.q
        delta_psi = (flux_in - q_seg / s.cell.material.xs * (1 - s.exponential[p]))
        return delta_psi

    def solveFlux(self, num_iter_tot, tol, update_source):

        stp = self.setup
        self.preCalculate()
        num_iter = 0
        delta_flux = 0
        #fuel = self.regions[0]
        #mod = self.regions[1]
        converged = False
        flux_old = []
        scalar_flux_old = [0, 0]

        print "Solving for fluxes...\n"

        # initialize scalar flux guesses for source calculation: phi = q / sigma_absorption for mod, phi = q for fuel
        if update_source:
            for i in range(self.mesh.n_cells):
                for cell in self.cells[i]:
                    #cell = self.cells[i][j]
                    if cell.material.name == 'fuel':
                        absorpt = 1
                    else:
                        absorpt = cell.material.xs - cell.material.scatter
                    cell.flux = cell.material.q / absorpt

        while not converged:
            if update_source:
                self.updateAllCellSource()
            for p in range(stp.n_p):                                       #loop over polar angles
                for i in range(stp.num_azim2):                             #loop over azimuthal angles
                    for track in self.tracks[i]:                           #loop over all tracks

                        track.quadwt = stp.omega_m[i] * stp.t_eff[i] * stp.omega_p[p] * stp.sintheta_p[p]  * 4 * math.pi
                        flux = track.flux_in[0, p]

                        if num_iter == 0:                                  #initial flux on boundaries is zero
                            flux = 0

                        #print "Forward tracking:"
                        for s in track.segments:                           #loop over segments
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux

                            #print "flux out %.1e \t\t delta flux %.1e" % (flux, delta_flux)

                            #increment flux in cell this segment crosses
                            s.cell.flux += delta_flux * track.quadwt / 2
                        #pass flux to next track
                        track.track_out.flux_in[track.refl_out, p] = flux

                        #track in reverse now
                        flux = track.flux_in[1,p]
                        if num_iter == 0:
                            flux = 0

                        #print "Reverse tracking:"
                        for s in track.segments[::-1]:
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux

                            # print " flux out %.1e \t\t delta flux %.1e" % (flux, delta_flux)
                            # increment flux in cell this segment crosses
                            s.cell.flux += delta_flux * track.quadwt / 2
                        # pass flux to next track
                        track.track_in.flux_in[track.refl_in, p] = flux

            #update fluxes in each FSR
            for k in range(self.mesh.n_cells):
                for cell in self.cells[k]:
                    cell.flux = (4 * math.pi * cell.material.q / cell.material.xs) \
                                + cell.flux / (cell.material.xs * cell.area)
                    #may need to have 4pi on the last term as well, check later

            getfluxes = self.getAvgScalarFlux()
            scalar_flux = getfluxes[:2]
            print "Checking convergence for iteration %d\n" % (num_iter)
            converged = self.check.isConverged(scalar_flux, scalar_flux_old, tol)
            #converged = True
            if not converged:
                num_iter += 1
                scalar_flux_old = scalar_flux[:]

                # reinitialize scalar flux accumulators in each FSR (cell)
                for k in range(self.mesh.n_cells):
                    for cell in self.cells[k]:
                        cell.flux = 0

                if num_iter == num_iter_tot:
                    converged = True
                    print "Not converged after %d iterations.\n" %(num_iter)
            else:
                print "Converged in %d iterations\n" %(num_iter)
                self.results = self.returnSolveResults(num_iter, getfluxes[0], getfluxes[1], getfluxes[2], getfluxes[3])


        #normalize fuel flux to 1
        #fuel.flux /= mod.flux
       # mod.flux /= mod.flux

        #mod.flux /= fuel.flux
        #fuel.flux /= fuel.flux

        #print "\nSCALAR FLUX\n-----------" \
        #      "\nFuel = \t\t\t%g \nModerator = \t%g" \
        #      "\nNumber of iterations: %d" % (fuel.flux, mod.flux, num_iter+1)

        #stp.plotScalarFlux(fuel.flux, mod.flux)
    def returnSolveResults(self, iters, fuelflux, modflux, avgflux, ratio):
       # print "Iterations to convergence: %d \nModerator flux: %g\nFuel flux: %g\n Avg flux: %g\n Flux ratio: %g" \
        #      %(iters, modflux, fuelflux, avgflux, ratio)
        return iters, fuelflux, modflux, avgflux, ratio

    def getAvgScalarFlux(self):
        fuelflux = 0.0
        modflux = 0.0
        maxflux = 0.0
        scalarflux = 0.0
        fuelcell = 0
        modcell = 0

        for i in range(self.mesh.n_cells):
            for cell in self.cells[i]:
                scalarflux += cell.flux
        print "scalar flux = %g" %(scalarflux)
        print "num cells = %g" %(self.mesh.n_cells ** 2)

        avg = scalarflux / (self.mesh.n_cells ** 2)
        for i in range(self.mesh.n_cells):
            for cell in self.cells[i]:

                if cell.flux > maxflux:
                    maxflux = cell.flux
                """
                if cell.region == 'fuel':
                    # accumulate scalar flux avg for fuel
                    fuelflux += cell.flux
                    fuelcell += 1
                else:
                    # accumulate scalar flux avg for mod
                    modflux += cell.flux
                    modcell += 1

                cell.flux /= avg
                """
        #avg_fuel = fuelflux / fuelcell
        #avg_mod = modflux / modcell


        for i in range(self.mesh.n_cells):
            for cell in self.cells[i]:
                cell.flux /= maxflux

                if cell.region == 'fuel':

                    #cell.flux /= avg_fuel

                    # accumulate scalar flux avg for fuel
                    fuelflux += cell.flux
                    fuelcell += 1

                elif cell.region == 'moderator':

                    #cell.flux /= avg_mod

                    # accumulate scalar flux avg for mod
                    modflux += cell.flux
                    modcell += 1

                else:
                    print "error in flux accumulation"

        fuelflux /= fuelcell
        print "max flux = %g" %(maxflux)
        #fuelflux /= maxflux
        modflux /= modcell
        print "fuelcell %g \t modcell %g\n" %(fuelcell, modcell)
        #modflux /= maxflux
        ratio = fuelflux / modflux
        avg /= maxflux
        print "Avg fuel flux = %f \nAvg mod flux = %f \nAverage Flux  = %f \nFlux ratio = %f" % (fuelflux, modflux, avg, ratio)
        return fuelflux, modflux, avg, ratio

    def updateSource(self, q, flux, scatter):
        #calc = ((1/2) * (scatter * flux + q))
        return ((1./4*math.pi) * (scatter * flux + q))

    def updateAllCellSource(self):
        for i in range(self.mesh.n_cells):
            for cell in self.cells[i]:
                #cell = self.cells[i][j]
               # print "cell [%d][%d] source update:\nOld source; %g" % (i, j, cell.source)

                cell.source = self.updateSource(cell.material.q, cell.flux, cell.material.scatter)
                #sourceval = cell.source
                #print "new source; %g" %(cell.source)

class ConvergenceTest(object):
    def __init__(self):
        """
        class for testing convergence using the L2 engineering norm
        n is the iteration index, i is the vector content index, I is total number of entries in vector.
        """

    def isConverged(self, vec_n, vec_n1, epsilon):
        sum1 = 0
        for i in range(len(vec_n)):
            error1 = ((vec_n[i] - vec_n1[i])/vec_n[i]) ** 2
            sum1 += error1
        I = len(vec_n)
        l2 = math.sqrt(sum1 / I)

        if l2 < epsilon:
            print "Converged! l2 %g" %(l2)
            return True
        else:
            print "Not converged; l2 %g" %(l2)
            return False

    def sourceProptoXSTest(self, xsfuel, xsmod):
        """test that sets source in each region proportional to the cross section in that region
        should yield a flat flux profile if code is functioning correctly"""
        qfuel = 2 * xsfuel
        qmod = 2 * xsmod
        print "This test case should yield a flat flux profile."
        return qfuel, qmod

    def sourceXSConstTest(self, qfuel, sigma_fuel):
        # angular flux everywhere should equal q/sigma
        qmod = qfuel
        sigma_mod = sigma_fuel
        print "Angular flux should equal %f everywhere when converged" %(qmod/sigma_mod)
        return qmod, sigma_mod

    def dancoffFactor(self, qfuel):
        qfuel = qfuel
        qmod = 0
        sigma_fuel = 1e5
        return qfuel, qmod, sigma_fuel

    def computeDancoff(self, phi_1, phi_fin, source, xs):
        const = 4 * math.pi * source / xs
        return 1 - (const - phi_fin)/ (const - phi_1)


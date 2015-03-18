#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import os
import numpy as np
import unittest

import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetect
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.utils.tests as utilsTests
import lsst.meas.algorithms.measureCurveOfGrowth as measComp

import lsst.afw.display.ds9 as ds9
try:
    display
except NameError:
    display = False

class CurveOfGrowthTestCase(utilsTests.TestCase):

    def setUp(self):
        self.config = measAlg.SourceMeasurementConfig()

        self.config.algorithms.names = ["flux.psf",
                                        "flux.aperture",
                                        "classification.extendedness"]
        self.config.algorithms["flux.aperture"].maxSincRadius = 10.0
        self.config.algorithms["flux.aperture"].radii = [3.0, 6.0, 9.0, 12.0, 15.0, 40, 65.0]
        self.config.slots.psfFlux = "flux.psf"

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField(afwTable.Field[int]("deblend.nchild", ""))

        task = measAlg.SourceMeasurementTask(config=self.config, schema=schema) # adds to schema

        sources = afwTable.SourceCatalog(schema)
        sources.definePsfFlux("flux.psf")
        sources.setMetadata(dafBase.PropertyList())
        sources.getMetadata().set("flux_aperture_radii", self.config.algorithms["flux.aperture"].radii)

        self.schema = schema
        self.sources = sources

    def tearDown(self):
        del self.schema
        del self.sources

    def fakeSources(self):
        """Add some fake sources to self.sources"""

        radii = list(self.config.algorithms["flux.aperture"].radii)
        area = np.pi*(np.array(radii + [0])**2 - np.array([0] + radii)**2)[0:-1]

        self.trueApertureFlux = np.cumsum(np.array([10, 5, 3, 1, 0.5, 0.1, 1e-10])*area)
        
        sky = 100
        for alpha in [1.0, 2.0, 3.0]:   # relative amplitudes of objects
            s = self.sources.addNew()
            s.set("flux.aperture", alpha*10*self.trueApertureFlux)
            s.set("flux.aperture.err", np.sqrt(sky*area + s.get("flux.aperture")))

            nInterpolatedPixel = np.zeros(len(radii), dtype='int32')
            if alpha == 1.0:
                nInterpolatedPixel[0] = 10 # one bad point
            s.set("flux.aperture.nInterpolatedPixel", nInterpolatedPixel)
                
            s.set("flux.aperture.nProfile", len(radii))
            s.set("flux.psf", alpha*1e8)

    def testCurveOfGrowth(self):
        """Test that CurveOfGrowthMeasurementTask generates a valid curve of growth
        """

        curveOfGrowthMeasurementTask = measComp.CurveOfGrowthMeasurementTask(None)

        self.fakeSources()
        result = curveOfGrowthMeasurementTask.run(self.sources)
        curveOfGrowth = result.curveOfGrowth
        #
        # Check that the curve of growth is the same as the input aperture fluxes
        #
        ratio = curveOfGrowth.apertureFlux/self.trueApertureFlux
        for r in ratio/ratio[0]:
            self.assertAlmostEqual(r, 1.0)

        if display:
            fig = curveOfGrowth.plot()
            if fig:
                fig.show()              # not in plot() for ipython's sake
                raw_input("Continue? ")
    
    def testPsfFluxFromCurveOfGrowth(self):
        """Test that we can estimate a psfFlux using the curveOfGrowth
        """

        curveOfGrowthMeasurementTask = measComp.CurveOfGrowthMeasurementTask(None)

        self.fakeSources()
        result = curveOfGrowthMeasurementTask.run(self.sources)
        curveOfGrowth = result.curveOfGrowth

        for s in self.sources:
            psfFlux = curveOfGrowth.estimatePsfFluxFromApertureFlux(s, i0=0, rMax=None)
            self.assertAlmostEqual(s.getPsfFlux()/psfFlux, 1)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CurveOfGrowthTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

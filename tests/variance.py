#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 LSST/AURA
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

"""
Tests for Variance measurement algorithm
"""

import unittest
import lsst.utils.tests as utilsTests

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as measAlg

try:
    display
except NameError:
    display = False


class VarianceTest(unittest.TestCase):
    def testVariance(self):
        size = 128 # size of image (pixels)
        center = afwGeom.Point2D(size//2, size//2) # object center
        width = 2.0 # PSF width
        flux = 10.0 # Flux of object
        variance = 1.0 # Variance value

        # Initial setup of an image
        exp = afwImage.ExposureF(size, size)
        image = exp.getMaskedImage().getImage()
        mask = exp.getMaskedImage().getMask()
        var = exp.getMaskedImage().getVariance()
        image.set(0.0)
        mask.set(0)
        var.set(variance)

        # Put down a PSF
        psf = measAlg.DoubleGaussianPsf(int(5*width), int(5*width), width)
        exp.setPsf(psf)
        psfImage = psf.computeImage(center).convertF()
        psfImage *= flux
        image.Factory(image, psfImage.getBBox(afwImage.PARENT)).__iadd__(psfImage)
        var.Factory(var, psfImage.getBBox(afwImage.PARENT)).__iadd__(psfImage)

        # Put in some bad pixels to ensure they're ignored
        for i in range(-5, 6):
            bad = size//2 + i*width
            var.getArray()[bad, :] = float("nan")
            mask.getArray()[bad, :] = mask.getPlaneBitMask("BAD")
            var.getArray()[:, bad] = float("nan")
            mask.getArray()[:, bad] = mask.getPlaneBitMask("BAD")

        if display:
            import lsst.afw.display.ds9 as ds9
            ds9.mtv(image, frame=1)
            ds9.mtv(mask, frame=2)
            ds9.mtv(var, frame=3)

        config = measAlg.SourceMeasurementConfig()
        config.algorithms.names = ["centroid.naive", "shape.sdss", "variance"]
        config.slots.centroid = "centroid.naive"
        config.slots.psfFlux = None
        config.slots.apFlux = None
        config.slots.modelFlux = None
        config.slots.instFlux = None
        config.slots.calibFlux = None
        config.slots.shape = "shape.sdss"
        config.algorithms["variance"].mask = ["BAD", "SAT"]

        config.validate()
        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = config.makeMeasureSources(schema)
        catalog = afwTable.SourceCatalog(schema)
        config.slots.setupTable(catalog.getTable())

        foot = afwDetection.Footprint(afwGeom.Point2I(center), width)
        peak = foot.getPeaks().addNew()
        peak.setIx(int(center.getX()))
        peak.setIy(int(center.getY()))
        peak.setFx(center.getX())
        peak.setFy(center.getY())
        peak.setPeakValue(flux)

        source = catalog.addNew()
        source.setFootprint(foot)
        ms.applyWithPeak(source, exp)

        self.assertEqual(source.get("variance"), variance)


##############################################################################################################

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(VarianceTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

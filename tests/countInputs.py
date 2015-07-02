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
Tests for CountInputs measurement algorithm
"""

import unittest
import lsst.utils.tests as utilsTests

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as measAlg
from lsst.afw.geom.polygon import Polygon


class CountInputsTest(unittest.TestCase):
    def testCountInputs(self):
        num = 3 # Number of images
        size = 10 # Size of images (pixels)
        shift = 4 # Shift to apply between images (pixels)
        value = 100.0 # Value to give objects
        offset = 12345 # x0,y0

        cdMatrix = (1.0e-5, 0.0, 0.0, 1.0e-5)
        crval = afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees)

        positions = [afwGeom.Point2D(size//2 + shift*(i - num//2) + offset,
                                     size//2 + shift*(i - num//2) + offset) for i in range(num)]
        wcsList = [afwImage.makeWcs(crval, pos - afwGeom.Extent2D(offset, offset), *cdMatrix) for
                   pos in positions]
        imageBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(size, size))
        wcsRef = afwImage.makeWcs(crval, afwGeom.Point2D(offset, offset), *cdMatrix)

        exp = afwImage.ExposureF(size, size)
        exp.setXY0(afwGeom.Point2I(offset, offset))
        exp.setWcs(wcsRef)
        exp.setPsf(measAlg.DoubleGaussianPsf(5, 5, 1.0))
        exp.getInfo().setCoaddInputs(afwImage.CoaddInputs(afwTable.ExposureTable.makeMinimalSchema(),
                                                          afwTable.ExposureTable.makeMinimalSchema()))
        ccds = exp.getInfo().getCoaddInputs().ccds
        for wcs in wcsList:
            record = ccds.addNew()
            record.setWcs(wcs)
            record.setBBox(imageBox)
            record.setValidPolygon(Polygon(afwGeom.Box2D(imageBox)))

        exp.getMaskedImage().getImage().set(0)
        exp.getMaskedImage().getMask().set(0)
        exp.getMaskedImage().getVariance().set(1.0)
        for pp in positions:
            x, y = map(int, pp)
            exp.getMaskedImage().getImage().set0(x, y, value)
            exp.getMaskedImage().getMask().set(x - offset, y - offset, value)

        measureSourcesConfig = measAlg.SourceMeasurementConfig()
        measureSourcesConfig.algorithms.names = ["centroid.naive", "countInputs"]
        measureSourcesConfig.slots.centroid = "centroid.naive"
        measureSourcesConfig.slots.psfFlux = None
        measureSourcesConfig.slots.apFlux = None
        measureSourcesConfig.slots.modelFlux = None
        measureSourcesConfig.slots.instFlux = None
        measureSourcesConfig.slots.calibFlux = None
        measureSourcesConfig.slots.shape = None
        measureSourcesConfig.validate()
        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = measureSourcesConfig.makeMeasureSources(schema)
        catalog = afwTable.SourceCatalog(schema)
        measureSourcesConfig.slots.setupTable(catalog.getTable())

        for pp in positions:
            foot = afwDetection.Footprint(afwGeom.Point2I(pp), 1.0)
            peak = foot.getPeaks().addNew()
            peak.setIx(int(pp[0]))
            peak.setIy(int(pp[1]))
            peak.setFx(pp[0])
            peak.setFy(pp[1])
            peak.setPeakValue(value)

            source = catalog.addNew()
            source.setFootprint(foot)

            ms.applyWithPeak(source, exp)

            number = sum(afwGeom.Box2D(imageBox).contains(wcs.skyToPixel(wcsRef.pixelToSky(pp))) for
                         wcs in wcsList)
            self.assertEqual(source.get("countInputs"), number)


##############################################################################################################

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CountInputsTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

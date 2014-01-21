#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

import numpy
import pickle
import unittest
import lsst.utils.tests as utilsTests

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.afw.geom import AffineTransform

from lsst.meas.algorithms import DiscreteBackground, CartesianPolygon, VectorPoint

DEBUG = False

class DiscreteBackgroundTest(unittest.TestCase):
    def setUp(self):
        self.x0, self.y0 = 1234, 5678             # Image offset from parent
        self.size = 250                           # Image size
        self.xCenter, self.yCenter = 123.4, 134.5 # Polygon center, in "LOCAL" frame
        self.xSize, self.ySize = 100, 200         # Polygon size
        self.offset = 1.234                       # Coefficient offset

    def polygon(self, angle):
        """Generate a polygon oriented at the nominated angle"""
        x0 = self.x0 + self.xCenter
        y0 = self.y0 + self.yCenter
        box = afwGeom.Box2D(afwGeom.Point2D(x0 - self.xSize/2, y0 - self.ySize/2),
                            afwGeom.Point2D(x0 + self.xSize/2, y0 + self.ySize/2))
        transform = (AffineTransform.makeTranslation(afwGeom.Extent2D(x0, y0)) *
                     AffineTransform.makeRotation(angle.asRadians())*
                     AffineTransform.makeTranslation(afwGeom.Extent2D(-x0, -y0)))
        poly = CartesianPolygon(box, transform)
        return poly

    def checkBackground(self, num):
        """Check background subtraction with 'num' polygons"""
        polyList = [self.polygon(afwGeom.Angle(afwGeom.PI/num*j)) for j in range(num)]
        image = afwImage.MaskedImageF(self.size, self.size)
        image.setXY0(self.x0, self.y0)
        image.set(0)
        coeffs = numpy.array([j + self.offset for j in range(num)])
        for c, p in zip(coeffs, polyList):
            polyImage = p.createImage(image.getBBox(afwImage.PARENT))
            image.getImage().scaledPlus(c, polyImage)
        bg = DiscreteBackground(image, polyList, 0)

        self.assertEqual(bg, bg)
        self.assertEqual(pickle.loads(pickle.dumps(bg)), bg)

        bgImage = bg.getImage()

        diff = image.clone()
        diff -= bgImage

        if DEBUG:
            import lsst.afw.display.ds9 as ds9
            ds9.mtv(image, frame=1, title="Image from %d polygons" % num)
            ds9.mtv(bgImage, frame=2, title="Background model with %d polygons" % num)
            ds9.mtv(diff, frame=3, title="Difference for %d polygons" % num)
            for frame in range(1, 4):
                for p in polyList:
                    p.display(frame=frame, xy0=image.getXY0())

        # Check coefficients
        solution = bg.getSolution()
        solRms = numpy.sqrt(numpy.average((solution - coeffs)**2))
        self.assertAlmostEqual(solRms, 0.0, 6)

        # Check subtracted image
        diffRms = numpy.sqrt(numpy.average(diff.getImage().getArray()**2))
        self.assertAlmostEqual(diffRms, 0.0, 5)

        print num, solRms, diffRms
        return bg

    def testBackground(self):
        bgList = []
        for num in range(3, 10):
            bg = self.checkBackground(num)
            for other in bgList:
                self.assertNotEqual(bg, other)
            bgList.append(bg)

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(DiscreteBackgroundTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

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
        self.x0 = 1234
        self.y0 = 5678
        self.xCenter = 123.4
        self.yCenter = 567.8
        self.xSize = 1000
        self.ySize = 2000
        self.offset = 1.234

    def polygon(self, angle):
        """Generate a polygon"""
        x0 = self.x0 + self.xCenter
        y0 = self.y0 + self.yCenter
        box = afwGeom.Box2D(afwGeom.Point2D(x0 - self.xSize/2, y0 - self.ySize/2),
                            afwGeom.Point2D(x0 + self.xSize/2, y0 + self.ySize/2))
        transform = (AffineTransform.makeTranslation(afwGeom.Extent2D(x0, y0)) *
                     AffineTransform.makeRotation(angle.asRadians())*
                     AffineTransform.makeTranslation(afwGeom.Extent2D(-x0, -y0)))
        return CartesianPolygon(box, transform)

    def testBackground(self):
        for i, num in enumerate(range(3, 10)):
            polyList = [self.polygon(afwGeom.Angle(2*afwGeom.PI/num*i)) for i in range(num)]
            image = afwImage.MaskedImageF(self.xSize, self.ySize)
            image.setXY0(self.x0, self.y0)
            image.set(0)
            for j, p in enumerate(polyList):
                polyImage = p.createImage(image.getBBox(afwImage.PARENT))
                image.getImage().scaledPlus(j + self.offset, polyImage)
            bg = DiscreteBackground(image, polyList, 0)
            bgImage = bg.getImage()

            solution = bg.getSolution()
            for j in range(num):
                # XXX may not be a good idea due to degeneracy
                self.assertAlmostEqual(solution[j], j + self.offset)

            diff = image.clone()
            diff -= bgImage

            if DEBUG:
                numPlots = 3
                import lsst.afw.display.ds9 as ds9
                ds9.mtv(image, frame=numPlots*i+1, title="Image from %d polygons" % num)
                ds9.mtv(bgImage, frame=numPlots*i+2, title="Background model with %d polygons" % num)
                ds9.mtv(diff, frame=numPlots*i+3, title="Difference for %d polygons" % num)
                for frame in range(numPlots*i+1, numPlots*(i+1)+1):
                    for p in polyList:
                        p.display(frame=frame)

            diff *= diff # chi^2
            print diff.sum()


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

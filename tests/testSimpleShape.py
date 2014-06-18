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

import unittest
import numpy
import os

import lsst.utils.tests
import lsst.afw.geom
import lsst.afw.geom.ellipses as el
import lsst.afw.image
from lsst.meas.algorithms import SimpleShape

class SimpleShapeTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.ellipseCores = [
            el.Quadrupole(25.0, 25.0, 0.0),
            el.Quadrupole(27.0, 22.0, -5.0),
            el.Quadrupole(23.0, 28.0, 2.0),
            ]
        self.centers = [
            lsst.afw.geom.Point2D(0.0, 0.0),
            lsst.afw.geom.Point2D(2.0, 3.0),
            lsst.afw.geom.Point2D(-1.0, 2.5),
            ]
        for ellipseCore in self.ellipseCores:
            ellipseCore.scale(2)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-500, -500), lsst.afw.geom.Point2I(50, 50))
        self.xg, self.yg = numpy.meshgrid(
            numpy.arange(self.bbox.getBeginX(), self.bbox.getEndX(), dtype=float),
            numpy.arange(self.bbox.getBeginY(), self.bbox.getEndY(), dtype=float)
            )

    def evaluateGaussian(self, ellipse):
        gt = ellipse.getGridTransform()
        xt = gt[gt.XX] * self.xg + gt[gt.XY] * self.yg + gt[gt.X]
        yt = gt[gt.YX] * self.xg + gt[gt.YY] * self.yg + gt[gt.Y]
        return numpy.exp(-0.5 * (xt**2 + yt**2))

    def checkMoments(self, dEllipseCore, dCenter, wEllipseCore, wCenter):
        dEllipse = el.Ellipse(dEllipseCore, dCenter)
        image = lsst.afw.image.MaskedImageF(self.bbox)
        image.getImage().getArray()[:,:] = self.evaluateGaussian(dEllipse)
        wEllipse = el.Ellipse(wEllipseCore, wCenter)
        result = SimpleShape.measure(wEllipse, image)
        return result

    def testCorrectWeightedMoments(self):
        for dEllipseCore in self.ellipseCores:
            for dCenter in self.centers:
                dEllipse = el.Ellipse(dEllipseCore, dCenter)
                dArray = self.evaluateGaussian(dEllipse)
                for wEllipseCore in self.ellipseCores:
                    for wCenter in self.centers:
                        wEllipse = el.Ellipse(wEllipseCore, wCenter)
                        wArray = self.evaluateGaussian(wEllipse)
                        product = dArray * wArray
                        i0 = numpy.sum(product)
                        ix = numpy.sum(product * self.xg) / i0
                        iy = numpy.sum(product * self.yg) / i0
                        ixx = numpy.sum(product * (self.xg - ix)**2) / i0
                        iyy = numpy.sum(product * (self.yg - iy)**2) / i0
                        ixy = numpy.sum(product * (self.xg - ix) * (self.yg - iy)) / i0
                        mEllipseCore = el.Quadrupole(ixx, iyy, ixy)
                        mCenter = lsst.afw.geom.Point2D(ix, iy)
                        SimpleShape.correctWeightedMoments(wEllipseCore, mEllipseCore, mCenter)
                        self.assertClose(mEllipseCore.getParameterVector(),
                                         dEllipseCore.getParameterVector(),
                                         rtol=1E-8, atol=1E-11)

    def testNoNoiseGaussians(self):
        for ellipseCore in self.ellipseCores:
            for center in self.centers:
                result = self.checkMoments(ellipseCore, center, ellipseCore, center)
                self.assertClose(result.ellipse.getParameterVector(),
                                 ellipseCore.getParameterVector(),
                                 rtol=3E-3, atol=1E-15)
                self.assertClose(numpy.array(result.center),
                                 numpy.array(center),
                                 rtol=1E-8, atol=1E-15)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SimpleShapeTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

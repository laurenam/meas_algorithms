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

    def polygon(self, angle):
        """Generate a polygon"""
        pass

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

#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2015 LSST/AURA
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
import unittest

import lsst.utils.tests
import lsst.meas.algorithms


class MomentsClassifierTestCase(lsst.utils.tests.TestCase):

    def testFullCoefficientArray(self):
        """Test that we can initialize coefficient arrays for specific filters as well as a backup."""
        for b in ((), ('r',), ('i',)):
            coefficients = numpy.random.randn(2, 3, 4)
            config = lsst.meas.algorithms.MomentsClassifierConfig()
            config.coefficients.clear()
            for i in xrange(coefficients.shape[0]):
                for j in xrange(coefficients.shape[1]):
                    for k in xrange(coefficients.shape[2]):
                        key = " ".join(b + tuple(str(v) for v in (i, j, k)))
                        config.coefficients[key] = float(coefficients[i, j, k])
            ctrl = config.makeControl()
            if b:
                self.assertClose(coefficients, ctrl.getCoefficients(b[0]), rtol=0.0, atol=0.0)
            else:
                self.assertClose(coefficients, ctrl.getCoefficients("q"), rtol=0.0, atol=0.0)

    def testSparseCoefficientArray(self):
        """Test that leaving out coefficients sets them to zero."""
        config = lsst.meas.algorithms.MomentsClassifierConfig()
        config.coefficients.clear()
        config.coefficients["0 0 2"] = 5.0
        ctrl = config.makeControl()
        self.assertClose(ctrl.getCoefficients("z"), numpy.array([[[0.0, 0.0, 5.0]]]), rtol=0.0, atol=0.0)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(MomentsClassifierTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests."""
    lsst.utils.tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

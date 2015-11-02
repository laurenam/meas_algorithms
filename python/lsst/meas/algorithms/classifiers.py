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

from . import algorithmsLib
from .algorithmRegistry import AlgorithmRegistry

def inClass(cls):
    """Decorator that allows a method to be defined as a free function."""
    def decorate(func):
        setattr(cls, func.__name__, func)
        return getattr(cls, func.__name__)
    return decorate


@inClass(algorithmsLib.MomentsClassifierConfig)
def setDefaults(self):
    super(algorithmsLib.MomentsClassifierConfig, self).setDefaults()
    # First index is trace of
    self.coefficients["0 0 0"] = -4.2759879274
    self.coefficients["0 1 0"] = 0.0713088756641
    self.coefficients["1 0 0"] = 0.16352932561
    self.coefficients["0 2 0"] = -4.54656639596e-05
    self.coefficients["1 1 0"] = -0.0482134274008
    self.coefficients["2 0 0"] = 4.41366874902e-13
    self.coefficients["0 3 0"] = 7.58973714641e-09
    self.coefficients["1 2 0"] = 1.51008430135e-05
    self.coefficients["2 1 0"] = 4.38493363998e-14
    self.coefficients["3 0 0"] = 1.83899834142e-20

AlgorithmRegistry.register(
    "classification.moments",
    target=algorithmsLib.MomentsClassifierControl,
    ConfigClass=algorithmsLib.MomentsClassifierConfig
)
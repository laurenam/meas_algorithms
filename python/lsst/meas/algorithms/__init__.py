# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

from .algorithmsLib import *
from .defects import *
from .psfDeterminerRegistry import *
from .starSelectorRegistry import *
from .algorithmRegistry import *
from .findCosmicRaysConfig import *
from .detection import *
from .measurement import *
from .replaceWithNoise import *
from .measureApCorr import *
from . import objectSizeStarSelector  # don't need names, just registration
from .makeCoaddApCorrMap import *
from .measureCurveOfGrowth import *

from .version import *

import lsst.utils

for name in dict(globals()):
    if name.endswith("_swigregister"):
        del globals()[name]

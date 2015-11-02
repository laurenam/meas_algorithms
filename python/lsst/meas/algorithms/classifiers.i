// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2015 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%{
#include "lsst/meas/algorithms/Extendedness.h"
#include "lsst/meas/algorithms/SigmoidClassifier.h"
#include "lsst/meas/algorithms/MomentsClassifier.h"
%}

%shared_ptr(lsst::meas::algorithms::ExtendednessControl)
%include "lsst/meas/algorithms/Extendedness.h"

%define %instantiateSigmoidClassifierPreInclude(N)
%shared_ptr(lsst::meas::algorithms::SigmoidClassifierControl<N>)
%shared_ptr(lsst::meas::algorithms::SigmoidClassifierAlgorithm<N>)
%declareNumPyConverters(ndarray::Array<double const,N,N>)
%enddef

%define %instantiateSigmoidClassifierPostInclude(N)
%template(SigmoidClassifierControl ## N) lsst::meas::algorithms::SigmoidClassifierControl<N>;
%template(SigmoidClassifierAlgorithm ## N) lsst::meas::algorithms::SigmoidClassifierAlgorithm<N>;
%enddef

%define %wrapSigmoidClassifier(NAME, N)
%pythoncode %{

def readSigmoidClassifierCoefficients(ctrl, flatDict):
    maxOrders = {}
    byFilter = {}
    for key, value in flatDict.iteritems():
        terms = key.split()
        if len(terms) == N:
            index = tuple(int(t) for t in terms)
            band = ""
        else:
            index = tuple(int(t) for t in terms[1:])
            band = terms[0]
        byFilter.setdefault(band, {})[index] = value
        for i in xrange(N):
            maxOrders[band, i] = max(maxOrders.get((band, i), 0), index[i])
    for band in byFilter:
        array = numpy.zeros(tuple((maxOrders[band, i] + 1) for i in xrange(N)), dtype=float)
        for index, value in byFilter[band].iteritems():
            array[index] = value
        ctrl.setCoefficients(band, array)

import lsst.pex.config
import numpy

@lsst.pex.config.wrap(NAME##Control)
class NAME##BaseConfig(AlgorithmConfig):

    coefficients = lsst.pex.config.DictField(
        keytype=str, itemtype=float, default={},
        doc=("Dictionary containing coefficients for the " #N "-d polynomial\n"
             "fed to the sigmoid function. Keys should be " #N " whitespace-\n"
             "delimited integers containing the index of the coefficient\n"
             "optionally preceded by the short form of the filter (i.e. 'i').\n"
             "Separate coefficient matrices may be defined for each filter, with\n"
             "an additional matrix used as the default.")
    )

class NAME ## Config(NAME##BaseConfig):

    def makeControl(self):
        ctrl = NAME##BaseConfig.makeControl(self)
        readSigmoidClassifierCoefficients(ctrl, self.coefficients)
        return ctrl

    def readControl(self, *args, **kwds):
        raise NotImplementedError("readControl() not implemented for SigmoidClassifiers")

%}
%enddef

%instantiateSigmoidClassifierPreInclude(1);
%instantiateSigmoidClassifierPreInclude(2);
%instantiateSigmoidClassifierPreInclude(3);

%include "lsst/meas/algorithms/SigmoidClassifier.h"

%instantiateSigmoidClassifierPostInclude(1);
%instantiateSigmoidClassifierPostInclude(2);
%instantiateSigmoidClassifierPostInclude(3);

%shared_ptr(lsst::meas::algorithms::MomentsClassifierControl)
%shared_ptr(lsst::meas::algorithms::MomentsClassifierAlgorithm)
%include "lsst/meas/algorithms/MomentsClassifier.h"

%wrapSigmoidClassifier(MomentsClassifier, 3)

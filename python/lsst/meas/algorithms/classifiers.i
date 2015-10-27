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
%enddef

%define %instantiateSigmoidClassifierPostInclude(N)
%template(SigmoidClassifierControl ## N) lsst::meas::algorithms::SigmoidClassifierControl<N>;
%template(SigmoidClassifierAlgorithm ## N) lsst::meas::algorithms::SigmoidClassifierAlgorithm<N>;
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
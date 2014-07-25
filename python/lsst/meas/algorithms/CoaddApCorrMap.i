// -*- lsst-C++ -*-

%{
#include "lsst/meas/algorithms/CoaddApCorrMap.h"
%}

%import "lsst/afw/image/ApCorrMap.i"

%shared_ptr(lsst::meas::algorithms::CoaddApCorrMap)
%declareTablePersistable(CoaddApCorrMap, lsst::meas::algorithms::CoaddApCorrMap);

%include "lsst/meas/algorithms/CoaddApCorrMap.h"

%castShared(lsst::meas::algorithms::CoaddApCorrMap, lsst::afw::image::ApCorrMap)

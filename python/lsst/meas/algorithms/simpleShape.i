// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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
#include "lsst/pex/logging.h"
#include "lsst/meas/algorithms.h"

#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_ALGORITHMS_SIMPLESHAPE_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%shared_ptr(lsst::meas::algorithms::SimpleShapeControl);
%shared_ptr(lsst::meas::algorithms::SimpleShape);

%include "lsst/meas/algorithms/SimpleShape.h"

%extend lsst::meas::algorithms::SimpleShape {
// add class-scope typedefs as typeobject class attributes, just
// to make C++ and Python interfaces more similar
%pythoncode %{
    Control = SimpleShapeControl
    Result = SimpleShapeResult
%}
}

%template(measure) lsst::meas::algorithms::SimpleShape::measure<float>;
%template(measure) lsst::meas::algorithms::SimpleShape::measure<double>;

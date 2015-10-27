// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2015 LSST/AURA
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

#include "boost/make_shared.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/meas/algorithms/MomentsClassifier.h"

namespace lsst { namespace meas { namespace algorithms {

MomentsClassifierControl::MomentsClassifierControl(std::string const & name) :
    SigmoidClassifierControl<3>(name),
    useDeterminant(false),
    useSdssBeforeHsm(true),
    flux("flux.psf")
{}

PTR(MomentsClassifierAlgorithm) MomentsClassifierControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    AlgorithmMap const & others,
    bool isForced
) const {
    return boost::static_pointer_cast<MomentsClassifierAlgorithm>(
        _makeAlgorithm(schema, metadata)
    );
}

PTR(AlgorithmControl) MomentsClassifierControl::_clone() const {
    return boost::make_shared<MomentsClassifierControl>(*this);
}

PTR(Algorithm) MomentsClassifierControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<MomentsClassifierAlgorithm>(*this, boost::ref(schema));
}

MomentsClassifierAlgorithm::MomentsClassifierAlgorithm(Control const & ctrl, afw::table::Schema & schema) :
    SigmoidClassifierAlgorithm<3>(
        ctrl,
        schema.addField<double>(
            ctrl.name,
            "star/galaxy classifier based a radius derived from a difference between source moments and "
            "PSF model moments."
        )
    ),
    _flux(schema[ctrl.flux]),
    _fluxErr(schema[ctrl.flux + ".err"]),
    _fluxFlag(schema[ctrl.flux + ".flags"]),
    _mainFlag(schema.addField<afw::table::Flag>(
        ctrl.name + ".flags", "Flag set when anything goes wrong with the classifier."
    )),
    _noShapeFlag(schema.addField<afw::table::Flag>(
        ctrl.name + ".flags.noShape", "Flag set if a needed shape input was unavailable or invalid."
    )),
    _noFluxFlag(schema.addField<afw::table::Flag>(
        ctrl.name + ".flags.noFlux", "Flag set if a needed flux input was unavailable or invalid."
    ))
{
    if (ctrl.useSdssBeforeHsm) {
        _shape1 = schema["shape.sdss"];
        _shapePsf1 = schema["shape.sdss.psf"];
        _shapeFlag1 = schema["shape.sdss.flags"];
        try {
            _shape2 = schema["shape.hsm.moments"];
            _shapePsf2 = schema["shape.hsm.psfMoments"];
            _shapeFlag2 = schema["shape.hsm.moments.flags"];
        } catch (pex::exceptions::NotFoundException &) {
            _shape2 = afw::table::Key< afw::table::Moments<double> >();
            _shapePsf2 = afw::table::Key< afw::table::Moments<double> >();
            _shapeFlag2 = afw::table::Key< afw::table::Flag >();
        }
    } else {
        _shape1 = schema["shape.hsm.moments"];
        _shapePsf1 = schema["shape.hsm.psfMoments"];
        _shapeFlag1 = schema["shape.hsm.moments.flags"];
        try {
            _shape2 = schema["shape.sdss"];
            _shapePsf2 = schema["shape.sdss.psf"];
            _shapeFlag2 = schema["shape.sdss.flags"];
        } catch (pex::exceptions::NotFoundException &) {
            _shape2 = afw::table::Key< afw::table::Moments<double> >();
            _shapePsf2 = afw::table::Key< afw::table::Moments<double> >();
            _shapeFlag2 = afw::table::Key< afw::table::Flag >();
        }
    }
}


void MomentsClassifierAlgorithm::getPoint(
    afw::table::SourceRecord & source,
    ndarray::Vector<double,3> & x,
    ndarray::Vector<bool,3> const & active,
    PTR(afw::image::Calib const) calib
) const {
    source.set(_mainFlag, true); // will unset later if we succeed
    double xx=0.0, yy=0.0, xy=0.0;
    if (!source.get(_shapeFlag1)) {
        xx = source.get(_shape1.getIxx()) - source.get(_shapePsf1.getIxx());
        yy = source.get(_shape1.getIyy()) - source.get(_shapePsf1.getIyy());
        xy = source.get(_shape1.getIxy()) - source.get(_shapePsf1.getIxy());
    } else if (_shapeFlag2.isValid() && !source.get(_shapeFlag2)) {
        xx = source.get(_shape2.getIxx()) - source.get(_shapePsf2.getIxx());
        yy = source.get(_shape2.getIyy()) - source.get(_shapePsf2.getIyy());
        xy = source.get(_shape2.getIxy()) - source.get(_shapePsf2.getIxy());
    } else {
        source.set(_noShapeFlag, true);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "No successful shape measurement available."
        );
    }
    if (getControl().useDeterminant) {
        x[0] = xx*yy - xy*xy;
    } else {
        x[0] = 0.5*(xx + yy);
    }
    if ((active[1] || active[2]) && source.get(_fluxFlag)) {
        source.set(_noFluxFlag, true);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "No successful flux measurement available."
        );
    }
    if (active[1]) {
        x[1] = source.get(_flux) / source.get(_fluxErr);
        if (!utils::isfinite(x[1])) {
            source.set(_noFluxFlag, true);
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Non-finite value in SNR."
            );
        }
    }
    if (active[2]) {
        source.set(_noFluxFlag, true); // in case we throw below (including in Calib.getMagnitude())
        if (!calib) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Calib required for magnitude-dependent classifier."
            );
        }
        x[2] = calib->getMagnitude(source.get(_flux));
        if (utils::isfinite(x[2])) {
            source.set(_noFluxFlag, false);
        } else {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Non-finite value in magnitude (probably a negative flux)."
            );
        }
    }
    source.set(_mainFlag, false); // unset the general failure flag we set earlier in case of exceptions
}

}}} // namespace lsst::meas::algorithms

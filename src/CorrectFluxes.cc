// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "lsst/utils/ieee.h"
#include "lsst/meas/algorithms/CorrectFluxes.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/image/ApCorrMap.h"

namespace lsst { namespace meas { namespace algorithms {

ApCorrRegistry & getApCorrRegistry() {
    static ApCorrRegistry instance;
    return instance;
}

namespace {

// Struct used to hold the Keys and field names associated with a particular flux field
// that needs to be aperture-corrected.
struct Element {

    Element(
        afw::table::Schema & schema,
        std::string const & fluxName_,
        bool doRecordApCorr
    ) : fluxName(fluxName_),
        fluxErrName(fluxName_ + ".err"),
        fluxKeys(schema[fluxName], schema[fluxErrName], schema[fluxName + ".flags"])
    {
        if (doRecordApCorr) {
            apCorrKeys.meas = schema.addField<afw::table::Flux::MeasTag>(
                fluxName + ".apcorr", "aperture correction applied to " + fluxName
            );
            apCorrKeys.err = schema.addField<afw::table::Flux::ErrTag>(
                fluxName + ".apcorr.err", "error on aperture correction applied to " + fluxName
            );
        }
        apCorrKeys.flag = schema.addField<afw::table::Flag>(
            fluxName + ".flags.apcorr", "set if aperture correction lookup failed"
        );
    }

    std::string fluxName;
    std::string fluxErrName;
    afw::table::KeyTuple<afw::table::Flux> fluxKeys;
    afw::table::KeyTuple<afw::table::Flux> apCorrKeys;
    mutable PTR(afw::math::BoundedField) apCorrField;
    mutable PTR(afw::math::BoundedField) apCorrErrField;
};

typedef std::vector<Element> ElementVector;

class CorrectFluxes : public Algorithm {
public:

    CorrectFluxes(
        CorrectFluxesControl const & ctrl,
        afw::table::Schema & schema,
        AlgorithmMap const & others
    );

private:

    void _apply(
        afw::table::SourceRecord & source,
        PTR(afw::image::ApCorrMap const) apCorrMap,
        afw::geom::Point2D const & center
    ) const;

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const {
        return this->_apply(source, exposure.getInfo()->getApCorrMap(), center);
    }

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(CorrectFluxes);

    ElementVector _elements;

    // We store the current ApCorrMap so we can check it against the one we're passed;
    // if they differ, we lookup all the fields and cache them so we don't have to look
    // them up for every source.
    mutable PTR(afw::image::ApCorrMap const) _apCorrMap;
};

CorrectFluxes::CorrectFluxes(
    CorrectFluxesControl const & ctrl,
    afw::table::Schema & schema,
    AlgorithmMap const & others
) : Algorithm(ctrl) {

    for (
        ApCorrRegistry::const_iterator nameIter = getApCorrRegistry().begin();
        nameIter != getApCorrRegistry().end();
        ++nameIter
    ) {
        if (std::find(ctrl.ignored.begin(), ctrl.ignored.end(), *nameIter) != ctrl.ignored.end()) {
            continue;
        }
        try {
            _elements.push_back(Element(schema, *nameIter, ctrl.doRecordApCorr));
        } catch (pex::exceptions::NotFoundException &) {
            // If the flux field field we want to correct isn't in the schema, we just ignore it.
        }
    }

}

void CorrectFluxes::_apply(
    afw::table::SourceRecord & source,
    PTR(afw::image::ApCorrMap const) apCorrMap,
    afw::geom::Point2D const & center
) const {
    CorrectFluxesControl const & ctrl = static_cast<CorrectFluxesControl const &>(this->getControl());

    if (!apCorrMap) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "No ApCorrMap attached to Exposure"
        );
    }

    if (apCorrMap != _apCorrMap) {
        // If we haven't set the cached ApCorrMap, or it doesn't match the one attached we've been passed,
        // we should update the per-(table-)field BoundedField objects.
        // This should happen once per exposure, on the first source processed.
        for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
            i->apCorrField = apCorrMap->get(i->fluxName);
            i->apCorrErrField = apCorrMap->get(i->fluxErrName);
        }
        _apCorrMap = apCorrMap;
    }

    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        // say we've failed when we start; we'll unset these flags when we succeed
        source.set(i->apCorrKeys.flag, true);
        bool oldFluxFlagState = false;
        if (ctrl.doFlagApCorrFailures) {
            oldFluxFlagState = source.get(i->fluxKeys.flag);
            source.set(i->fluxKeys.flag, true);
        }
        if (!i->apCorrField || !i->apCorrErrField) {
            continue;
        }
        double apCorr = 1.0;
        double apCorrErr = 0.0;
        try {
            apCorr = i->apCorrField->evaluate(center);
            if (ctrl.doRecordApCorr) {
                source.set(i->apCorrKeys.meas, apCorr);
            }
        } catch (pex::exceptions::Exception &) {
            continue;
        }
        try {
            // We attempt to record the aperture correction error, even if we don't propagate it into the
            // flux errors - but if we didn't measure the aperture correction error, it's only a fatal error
            // if we are trying to propagate it.
            apCorrErr = i->apCorrErrField->evaluate(center);
            if (ctrl.doRecordApCorr) {
                source.set(i->apCorrKeys.err, apCorrErr);
            }
        } catch (pex::exceptions::Exception &) {
            if (ctrl.doPropagateErrors) continue;
        }
        if (apCorr <= 0.0 || (ctrl.doPropagateErrors && apCorrErr < 0.0)) {
            continue;
        }
        double const flux = source.get(i->fluxKeys.meas);
        source.set(i->fluxKeys.meas, flux*apCorr);
        double const fluxErr = source.get(i->fluxKeys.err);
        if (ctrl.doPropagateErrors) {
            double const a = fluxErr/flux;
            double const b = apCorrErr/apCorr;
            source.set(i->fluxKeys.err, std::abs(flux*apCorr)*std::sqrt(a*a + b*b));
        } else {
            source.set(i->fluxKeys.err, apCorr*fluxErr);
        }
        source.set(i->apCorrKeys.flag, false);
        if (ctrl.doFlagApCorrFailures) {
            source.set(i->fluxKeys.flag, oldFluxFlagState);
        }
    }
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(CorrectFluxes);

} // anonymous

PTR(AlgorithmControl) CorrectFluxesControl::_clone() const {
    return boost::make_shared<CorrectFluxesControl>(*this);
}

PTR(Algorithm) CorrectFluxesControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    AlgorithmMap const & others,
    bool isForced
) const {
    return boost::make_shared<CorrectFluxes>(*this, boost::ref(schema), others);
}

}}} // namespace lsst::meas::algorithms

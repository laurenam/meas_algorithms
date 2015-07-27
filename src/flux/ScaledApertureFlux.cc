// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/FluxControl.h"

namespace lsst {
namespace meas {
namespace algorithms {

class ScaledApertureFlux : public FluxAlgorithm {
public:

    ScaledApertureFlux(ScaledApertureFluxControl const & ctrl, afw::table::Schema & schema) :
        FluxAlgorithm(ctrl, schema, "photometry with aperture scaled by PSF") {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(ScaledApertureFlux);

};

template <typename PixelT>
void ScaledApertureFlux::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    ScaledApertureFluxControl const& ctrl = static_cast<ScaledApertureFluxControl const&>(getControl());

    double const radius = exposure.getPsf()->computeShape(center).getDeterminantRadius();
    double const fwhm = 2.0*std::sqrt(2.0*std::log(2))*radius;
    double const size = ctrl.scale*fwhm;
    afw::geom::ellipses::Axes const axes(size, size);

    std::pair<double, double> fluxes =
        photometry::calculateSincApertureFlux(exposure.getMaskedImage(),
                                              afw::geom::ellipses::Ellipse(axes, center));
    source.set(getKeys().meas, fluxes.first);
    source.set(getKeys().err, fluxes.second);
    source.set(getKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(ScaledApertureFlux);

PTR(AlgorithmControl) ScaledApertureFluxControl::_clone() const {
    return boost::make_shared<ScaledApertureFluxControl>(*this);
}

PTR(Algorithm) ScaledApertureFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const&
) const {
    return boost::make_shared<ScaledApertureFlux>(*this, boost::ref(schema));
}


}}}

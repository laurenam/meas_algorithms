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

#ifndef LSST_MEAS_ALGORITHMS_Blendedness_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_Blendedness_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Image.h"

namespace lsst { namespace meas { namespace algorithms {


struct BlendednessControl {

    LSST_CONTROL_FIELD(
        doOld, bool,
        "Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)"
    );

    LSST_CONTROL_FIELD(
        doFlux, bool,
        "Whether to compute quantities related to the Gaussian-weighted flux"
    );

    LSST_CONTROL_FIELD(
        doShape, bool,
        "Whether to compute quantities related to the Gaussian-weighted shape"
    );

    LSST_CONTROL_FIELD(
        nSigmaWeightMax, double,
        "Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)"
    );

    LSST_CONTROL_FIELD(
        useAbsoluteValue, bool,
        "If true, take the absolute value of the pixels when computing flux and shape blendedness."
    );

    BlendednessControl() :
        doOld(true),
        doFlux(true),
        doShape(true),
        nSigmaWeightMax(3.0),
        useAbsoluteValue(true)
    {}

};

/**
 *  Compute metrics that measure how blended objects are.
 *
 *  Blendedness is initialized once for a given Schema, then run repeatedly
 *  by calls to measureChildPixels and measureParentPixels (in any order, possibly
 *  with multiple sources interleaved), followed by a call to finish() for each
 *  source.
 */
class Blendedness {
public:

    explicit Blendedness(BlendednessControl const & ctrl, afw::table::Schema & schema);

    void measureChildPixels(
        afw::image::Image<float> const & image,
        afw::table::SourceRecord & child
    ) const;

    void measureParentPixels(
        afw::image::Image<float> const & image,
        afw::table::SourceRecord & child
    ) const;

private:

    void _measureMoments(
        afw::image::Image<float> const & image,
        afw::table::SourceRecord & child,
        afw::table::Key<double> const & fluxKey,
        afw::table::Key< afw::table::Moments<double> > const & _shapeKey
    ) const;

    BlendednessControl const _ctrl;
    afw::table::Key<double> _old;
    afw::table::Key<double> _flux;
    afw::table::Key<double> _fluxChild;
    afw::table::Key<double> _fluxParent;
    afw::table::Key< afw::table::Moments<double> > _shapeChild;
    afw::table::Key< afw::table::Moments<double> > _shapeParent;
    afw::table::Key<afw::table::Flag> _flagGeneral;
    afw::table::Key<afw::table::Flag> _flagNoCentroid;
    afw::table::Key<afw::table::Flag> _flagNoShape;
};

}}} // namespace lsst::meas::algorithms

#endif // LSST_MEAS_ALGORITHMS_Blendedness_h_INCLUDED

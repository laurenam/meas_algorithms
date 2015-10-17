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

    BlendednessControl() :
        doOld(true)
    {}

};

/**
 *  Compute metrics that measure how blended objects are.
 *
 *  Blendedness is initialized once for a given Schema, then run repeatedly
 *  by calls to measureChildPixels and measParentPixels (in any order, possibly
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
    BlendednessControl const _ctrl;
    afw::table::Key<double> _old;
};

}}} // namespace lsst::meas::algorithms

#endif // LSST_MEAS_ALGORITHMS_Blendedness_h_INCLUDED

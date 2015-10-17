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

#include "lsst/meas/algorithms/Blendedness.h"
#include "lsst/afw/detection/HeavyFootprint.h"


namespace lsst { namespace meas { namespace algorithms {


namespace {


double computeOldBlendedness(
    PTR(afw::detection::Footprint const) childFootprint,
    afw::image::Image<float> const & parentImage
) {
    if (!childFootprint) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "blendedness.old requires a Footprint."
        );
    }

    PTR(afw::detection::HeavyFootprint<float> const) childHeavy =
        boost::dynamic_pointer_cast<afw::detection::HeavyFootprint<float> const>(childFootprint);

    if (!childHeavy) {
        return 0.0;  // if it's not a HeavyFootprint, it's not blended.
    }

    if (!parentImage.getBBox(afw::image::PARENT).contains(childHeavy->getBBox())) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Child footprint extends beyond image."
        );
    }

    // Iterate over all the spans in the child HeavyFootprint,
    // along with iterators for the child pixels (from the HeavyFootprint)
    // and parent pixels (from the Exposure).
    typedef afw::detection::Footprint::SpanList::const_iterator SpanIter;
    typedef afw::image::Image<float>::const_x_iterator ParentPixIter;
    typedef ndarray::Array<float const,1,1>::Iterator ChildPixIter;
    SpanIter spanIter = childHeavy->getSpans().begin();
    SpanIter const spanEnd = childHeavy->getSpans().end();
    ChildPixIter childPixIter = childHeavy->getImageArray().begin();
    double cp = 0.0;   // child.dot(parent)
    double cc = 0.0;   // child.dot(child)
    while (spanIter != spanEnd) {
        afw::geom::Span const & span = **spanIter;
        ParentPixIter parentPixIter = parentImage.x_at(
            span.getBeginX() - parentImage.getX0(),
            span.getY() - parentImage.getY0()
        );
        int const width = span.getWidth();
        // Iterate over the pixels within the span, updating the dot products.
        for (int x = 0; x < width; ++parentPixIter, ++childPixIter, ++x) {
            cp += (*childPixIter) * ((*parentPixIter) - (*childPixIter));
            cc += (*childPixIter) * (*childPixIter);
        }
        ++spanIter;
    }
    if (cc > 0.0) {
        return cp/cc;
    }
    return 0.0;
}

} // anonymous


Blendedness::Blendedness(BlendednessControl const & ctrl, afw::table::Schema & schema) :
    _ctrl(ctrl)
{
    if (_ctrl.doOld) {
        _old = schema.addField<double>(
            "blendedness.old",
            "blendedness from dot products: (child.dot(parent)/child.dot(child) - 1)"
        );
    }
}


void Blendedness::measureChildPixels(
    afw::image::Image<float> const & image,
    afw::table::SourceRecord & child
) const {
    // Nothing to do here for old blendedness
}


void Blendedness::measureParentPixels(
    afw::image::Image<float> const & image,
    afw::table::SourceRecord & child
) const {
    if (_ctrl.doOld) {
        child.set(_old, computeOldBlendedness(child.getFootprint(), image));
    }
}

}}} // namespace lsst::meas::algorithms

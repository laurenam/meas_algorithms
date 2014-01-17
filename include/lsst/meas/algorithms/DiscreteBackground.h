// -*- LSST-C++ -*-

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

#if !defined(LSST_MEAS_ALGORITHMS_DISCRETE_BACKGROUND_H)
#define LSST_MEAS_ALGORITHMS_DISCRETE_BACKGROUND_H

#include <vector>

#include "ndarray_fwd.h"

#include "lsst/base.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/math/Background.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"

#include "lsst/meas/algorithms/CartesianPolygon.h"

namespace lsst { namespace meas { namespace algorithms {

class DiscreteBackground {
public:
    typedef float PixelT;
    typedef afw::image::Exposure<PixelT, afw::image::MaskPixel> Exposure;
    typedef std::vector<afw::geom::Box2I> BoxVector;
    typedef std::vector<CONST_PTR(afw::image::Wcs)> WcsVector;
    typedef std::vector<CartesianPolygon> PolygonVector;
    typedef ndarray::Array<double const, 1, 1> Solution;

    DiscreteBackground(Exposure const& exp,
                       afw::geom::Box2I const& box,
                       CONST_PTR(afw::image::Wcs) const& wcs,
                       afw::image::MaskPixel const maskVal
        );
    DiscreteBackground(Exposure const& exp,
                       BoxVector const& boxList,
                       WcsVector const& wcsList,
                       afw::image::MaskPixel const maskVal
        );
    DiscreteBackground(PolygonVector const& polyList,
                       Solution const& solution,
                       afw::geom::Box2I const& bbox
        );
    DiscreteBackground(DiscreteBackground const& other) : _impl(other._impl) {}
    virtual ~DiscreteBackground() {}

    PTR(afw::image::Image<PixelT>) getImage() const;

    PolygonVector getPolygonVector() const;
    Solution getSolution() const;
    afw::geom::Box2I getBBox() const;

private:
    DiscreteBackground& operator=(CartesianPolygon const&); // assignment unimplemented

    /// pImpl pattern to hide implementation
    struct Impl;
    PTR(Impl) _impl;
    DiscreteBackground(PTR(Impl) impl) : _impl(impl) {}
};

}}} // namespace lsst::meas::algorithms

#endif

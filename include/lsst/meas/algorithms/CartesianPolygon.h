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

#if !defined(LSST_MEAS_ALGORITHMS_CARTESIAN_POLYGON_H)
#define LSST_MEAS_ALGORITHMS_CARTESIAN_POLYGON_H

#include <vector>
#include <utility> // for std::pair

#include "boost/make_shared.hpp"

#include "lsst/base.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image/Wcs.h"

namespace lsst { namespace meas { namespace algorithms {

class CartesianPolygon {
public:
    typedef afw::geom::Box2D Box;
    typedef afw::geom::Point2D Point;

    /// Constructors
    explicit CartesianPolygon(Box const& box);
    explicit CartesianPolygon(
        Box const &box,                  ///< Box to convert to polygon
        afw::image::Wcs const& original, ///< Original frame Wcs for Box
        afw::image::Wcs const& target    ///< Target frame Wcs for polygon
        );
    explicit CartesianPolygon(std::vector<Point> const& vertices);
    virtual ~CartesianPolygon() {}

    // Copying just involves copying the implementation
    CartesianPolygon(CartesianPolygon const& other) : _impl(other._impl) {}
    CartesianPolygon& operator=(CartesianPolygon const& other) {
        _impl = other._impl;
        return *this;
    }

    /// Return number of edges
    ///
    /// Identical with the number of points
    size_t getNumEdges() const;

    /// Return bounding box
    Box getBBox() const;

    Point calculateCenter() const;
    double calculateArea() const;
    double calculatePerimeter() const;

    /// Get vector of vertices
    ///
    /// Note that the "closed" polygon vertices are returned, so the first and
    /// last vertex are identical and there is one more vertex than otherwise
    /// expected.
    std::vector<Point> getVertices() const;

    /// Get vector of edges
    ///
    /// Returns edges, as pairs of vertices.
    std::vector<std::pair<Point, Point> > getEdges() const;

    bool operator==(CartesianPolygon const& other) const;
    bool operator!=(CartesianPolygon const& other) const { return !(*this == other); }

    /// Returns whether the polygon contains the point
    bool contains(Point const& point) const;

    /// Returns whether the polygons overlap each other
    bool overlaps(CartesianPolygon const& other) const;

    /// Returns the intersection of two polygons
    ///
    /// Does not handle non-convex polygons (which might have multiple independent
    /// intersections).
    CartesianPolygon intersection(CartesianPolygon const& other) const;

    CartesianPolygon convexHull() const;

    /// Create image of polygon
    ///
    /// Pixels entirely contained within the polygon receive value unity,
    /// pixels entirely outside the polygon receive value zero, and pixels
    /// on the border receive a value equal to the fraction of the pixel
    /// within the polygon.
    ///
    /// Note that the center of the lower-left pixel is 0,0.
    PTR(afw::image::Image<float>) createImage(afw::geom::Box2I const& bbox) const;
    PTR(afw::image::Image<float>) createImage(afw::geom::Extent2I const& extent) const {
        return createImage(afw::geom::Box2I(afw::geom::Point2I(0, 0), extent));
    }

private:
    /// pImpl pattern to hide implementation
    struct Impl;
    PTR(Impl) _impl;
    CartesianPolygon(PTR(Impl) impl) : _impl(impl) {}
};

std::ostream& operator<<(std::ostream& os, CartesianPolygon const& box);

}}} // namespace lsst::meas::algorithms

#endif

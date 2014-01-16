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
    CartesianPolygon(Box const& box);
    CartesianPolygon(
        Box const &box,                  ///< Box to convert to polygon
        afw::image::Wcs const& original, ///< Original frame Wcs for Box
        afw::image::Wcs const& target    ///< Target frame Wcs for polygon
        );
    CartesianPolygon(std::vector<Point> const& vertices);
    virtual ~CartesianPolygon() {}
    CartesianPolygon(CartesianPolygon const& other) : _impl(other._impl) {}

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

private:
    CartesianPolygon& operator=(CartesianPolygon const&); // assignment unimplemented

    /// pImpl pattern to hide implementation
    struct Impl;
    PTR(Impl) _impl;
    CartesianPolygon(PTR(Impl) impl) : _impl(impl) {}
};


}}} // namespace lsst::meas::algorithms

#endif

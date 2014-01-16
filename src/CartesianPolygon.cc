#include "boost/geometry/geometry.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/CartesianPolygon.h"



typedef lsst::meas::algorithms::CartesianPolygon::Point LsstPoint;
typedef lsst::meas::algorithms::CartesianPolygon::Box LsstBox;
typedef std::vector<LsstPoint> LsstRing;
typedef boost::geometry::model::polygon<LsstPoint> BoostPolygon;
typedef boost::geometry::model::box<LsstPoint> BoostBox;

namespace boost { namespace geometry { namespace traits {

// Setting up LsstPoint
template<> struct tag<LsstPoint> { typedef point_tag type; };
template<> struct coordinate_type<LsstPoint> { typedef LsstPoint::Element type; };
template<> struct coordinate_system<LsstPoint> { typedef cs::cartesian type; };
template<> struct dimension<LsstPoint> : boost::mpl::int_<2> {};
template<unsigned long dim>
struct access<LsstPoint, dim>
{
    static double get(LsstPoint const& p) { return p[dim]; }
    static void set(LsstPoint& p, LsstPoint::Element const& value) { p[dim] = value; }
};

// Setting up LsstBox
//
// No setters, because it's inefficient (can't set individual elements of Box2D directly).
// For box outputs from boost::geometry we'll use BoostBox and then convert.
template<> struct tag<LsstBox> { typedef box_tag type; };
template<> struct point_type<LsstBox> { typedef LsstPoint type; };
template<>
struct indexed_access<LsstBox, 0, 0>
{
    static double get(LsstBox const& box) { return box.getMinX(); }
};
template<>
struct indexed_access<LsstBox, 1, 0>
{
    static double get(LsstBox const& box) { return box.getMaxX(); }
};
template<>
struct indexed_access<LsstBox, 0, 1>
{
    static double get(LsstBox const& box) { return box.getMinY(); }
};
template<>
struct indexed_access<LsstBox, 1, 1>
{
    static double get(LsstBox const& box) { return box.getMaxY(); }
};

// Setting up LsstRing
template<> struct tag<LsstRing> { typedef ring_tag type; };
//template<> struct range_value<LsstRing> { typedef LsstPoint type; };

}}} // namespace boost::geometry::traits


namespace {

/// Convert BoostBox to LsstBox
LsstBox boostBoxToLsst(BoostBox const& box)
{
    return LsstBox(box.min_corner(), box.max_corner());
}

} // anonymous namespace


namespace lsst { namespace meas { namespace algorithms {

struct CartesianPolygon::Impl
{
    Impl() : poly() {}
    explicit Impl(CartesianPolygon::Box const& box) : poly() {
        boost::geometry::assign(poly, box);
        check();
    }
    explicit Impl(std::vector<CartesianPolygon::Point> const& vertices) : poly() {
        boost::geometry::assign(poly, vertices);
        check();
    }
    explicit Impl(BoostPolygon _poly) : poly(_poly) {}

    void check() { boost::geometry::correct(poly); }

    BoostPolygon poly;
};


CartesianPolygon::CartesianPolygon(CartesianPolygon::Box const& box) :
    _impl(new CartesianPolygon::Impl(box)) {}

CartesianPolygon::CartesianPolygon(std::vector<CartesianPolygon::Point> const& vertices) :
    _impl(new CartesianPolygon::Impl(vertices)) {}

CartesianPolygon::CartesianPolygon(
    CartesianPolygon::Box const &box,
    afw::image::Wcs const& original,
    afw::image::Wcs const& target
    ) : _impl(new CartesianPolygon::Impl())
{
    std::vector<lsst::afw::geom::Point2D> vertices;
    vertices.reserve(4);
    vertices.push_back(box.getMin());
    vertices.push_back(lsst::afw::geom::Point2D(box.getMaxX(), box.getMinY()));
    vertices.push_back(box.getMax());
    vertices.push_back(lsst::afw::geom::Point2D(box.getMinX(), box.getMaxY()));
    for (typename std::vector<lsst::afw::geom::Point2D>::iterator p = vertices.begin();
         p != vertices.end(); ++p) {
        *p = target.skyToPixel(*original.pixelToSky(*p));
    }
    boost::geometry::assign(_impl->poly, vertices);
    _impl->check();
}

size_t CartesianPolygon::getNumEdges() const {
    // boost::geometry::models::polygon uses a "closed" polygon: the start/end point is included twice
    return boost::geometry::num_points(_impl->poly) - 1;
}

CartesianPolygon::Box CartesianPolygon::getBBox() const
{
    return boostBoxToLsst(boost::geometry::return_envelope<BoostBox>(_impl->poly));
}

CartesianPolygon::Point CartesianPolygon::calculateCenter() const
{
    return boost::geometry::return_centroid<LsstPoint>(_impl->poly);
}

double CartesianPolygon::calculateArea() const { return boost::geometry::area(_impl->poly); }

double CartesianPolygon::calculatePerimeter() const { return boost::geometry::perimeter(_impl->poly); }

std::vector<std::pair<CartesianPolygon::Point, CartesianPolygon::Point> >
CartesianPolygon::getEdges() const
{
    std::vector<LsstPoint> const vertices = getVertices();
    std::vector<std::pair<CartesianPolygon::Point, CartesianPolygon::Point> > edges;
    edges.reserve(getNumEdges());
    for (typename std::vector<LsstPoint>::const_iterator i = vertices.begin(), j = vertices.begin() + 1;
         j != vertices.end(); ++i, ++j) {
        edges.push_back(std::make_pair(*i, *j));
    }
    return edges;
}

std::vector<CartesianPolygon::Point> CartesianPolygon::getVertices() const { return _impl->poly.outer(); }

bool CartesianPolygon::operator==(CartesianPolygon const& other) const
{
    return boost::geometry::equals(_impl->poly, other._impl->poly);
}

bool CartesianPolygon::contains(CartesianPolygon::Point const& point) const
{
    return boost::geometry::within(point, _impl->poly);
}

bool CartesianPolygon::overlaps(CartesianPolygon const& other) const {
    return !boost::geometry::disjoint(_impl->poly, other._impl->poly);
}

CartesianPolygon CartesianPolygon::intersection(CartesianPolygon const& other) const
{
    std::vector<BoostPolygon> result;
    boost::geometry::intersection(_impl->poly, other._impl->poly, result);
    if (result.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Polygons have no intersection");
    }
    if (result.size() > 1) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException,
                          (boost::format("Multiple polygons (%d) created by intersection()") %
                           result.size()).str());
    }
    return CartesianPolygon(PTR(Impl)(new Impl(result[0])));
}


}}} // namespace lsst::meas::algorithms

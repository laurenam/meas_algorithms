#include <cmath>
#include <algorithm>

#include "boost/geometry/geometry.hpp"
#include "boost/make_shared.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/CartesianPolygon.h"



typedef lsst::meas::algorithms::CartesianPolygon::Point LsstPoint;
typedef lsst::meas::algorithms::CartesianPolygon::Box LsstBox;
typedef std::vector<LsstPoint> LsstRing;
typedef boost::geometry::model::polygon<LsstPoint> BoostPolygon;
typedef boost::geometry::model::box<LsstPoint> BoostBox;
typedef boost::geometry::model::linestring<LsstPoint> BoostLineString;

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

/// Convert box to corners
std::vector<LsstPoint> boxToCorners(LsstBox const& box)
{
    std::vector<LsstPoint> corners;
    corners.reserve(4);
    corners.push_back(box.getMin());
    corners.push_back(LsstPoint(box.getMaxX(), box.getMinY()));
    corners.push_back(box.getMax());
    corners.push_back(LsstPoint(box.getMinX(), box.getMaxY()));
    return corners;
}

/// Calculate area of overlap between polygon and pixel
double pixelOverlap(BoostPolygon const& poly, int const x, int const y)
{
    BoostPolygon pixel;     // Polygon for single pixel
    boost::geometry::assign(pixel, LsstBox(lsst::afw::geom::Point2D(x - 0.5, y - 0.5),
                                           lsst::afw::geom::Extent2D(1.0, 1.0)));
    std::vector<BoostPolygon> overlap; // Overlap between pixel and polygon
    boost::geometry::intersection(poly, pixel, overlap);
    double area = 0.0;
    for (std::vector<BoostPolygon>::const_iterator i = overlap.begin(); i != overlap.end(); ++i) {
        area += boost::geometry::area(*i);
        assert(area <= 1.0); // by construction
    }
    return area;
}

/// Set each pixel in a row to the amount of overlap with polygon
void pixelRowOverlap(PTR(lsst::afw::image::Image<float>) const image,
                     BoostPolygon const& poly, int const xStart, int const xStop, int const y)
{
    int x = xStart;
    for (lsst::afw::image::Image<float>::x_iterator i = image->x_at(x - image->getX0(), y - image->getY0());
         x <= xStop; ++i, ++x) {
        *i = pixelOverlap(poly, x, y);
    }
}

template <typename PixelT>
struct DummyIterator {
    DummyIterator(int _value) : value(_value) {}
    PixelT operator*() const { return value; }
    DummyIterator const& operator++() const { return *this; } // prefix increment
    DummyIterator const& operator++(int) const { return *this; } // postfix increment
    PixelT const value;
};

// Duck-typing an Image with value identically unity
struct DummyImage {
    typedef DummyIterator<float> x_iterator;
    x_iterator x_at(int x, int y) const { return x_iterator(1.0); }
};

// Duck-typing a Mask with value identically zero
struct DummyMask {
    typedef DummyIterator<lsst::afw::image::MaskPixel> x_iterator;
    x_iterator x_at(int x, int y) const { return x_iterator(0); }
};

// Duck-typing a MaskedImage
template <typename ImageT, typename MaskT>
struct DummyMaskedImage {
    typedef MaskT Mask;
    typedef ImageT Image;
    DummyMaskedImage(CONST_PTR(ImageT) _image, CONST_PTR(MaskT) _mask, lsst::afw::geom::Box2I const& _bbox)
        : image(_image), mask(_mask), bbox(_bbox) {}
    CONST_PTR(Mask) getMask() const { return mask; }
    CONST_PTR(Image) getImage() const { return image; }
    lsst::afw::geom::Box2I const& getBBox(lsst::afw::image::ImageOrigin) const { return bbox; }
    int getX0() const { return bbox.getMin().getX(); }
    int getY0() const { return bbox.getMin().getY(); }
    CONST_PTR(ImageT) const image;
    CONST_PTR(MaskT) const mask;
    lsst::afw::geom::Box2I const& bbox;
};


} // anonymous namespace


namespace lsst { namespace meas { namespace algorithms {

/// Stream vertices
std::ostream& operator<<(std::ostream& os, std::vector<LsstPoint> const& vertices)
{
    os << "[";
    size_t num = vertices.size();
    for (size_t i = 0; i < num - 1; ++i) {
        os << vertices[i] << ",";
    }
    os << vertices[vertices.size() - 1] << "]";
    return os;
}

/// Stream BoostPolygon
std::ostream& operator<<(std::ostream& os, BoostPolygon const& poly)
{
    return os << "BoostPolygon(" << poly.outer() << ")";
}

std::ostream& operator<<(std::ostream& os, CartesianPolygon const& poly)
{
    os << "CartesianPolygon(" << poly.getVertices() << ")";
    return os;
}

struct CartesianPolygon::Impl
{
    Impl() : poly() {}
    explicit Impl(CartesianPolygon::Box const& box) : poly() {
        boost::geometry::assign(poly, box);
        // Assignment from a box is correctly handled by BoostPolygon, so doesn't need a "check()"
    }
    explicit Impl(std::vector<CartesianPolygon::Point> const& vertices) : poly() {
        boost::geometry::assign(poly, vertices);
        check(); // because the vertices might not have the correct orientation (CW vs CCW) or be open
    }
    explicit Impl(BoostPolygon const& _poly) : poly(_poly) {}

    void check() { boost::geometry::correct(poly); }

    template <class PolyT>
    bool overlaps(PolyT const& other) const {
        return !boost::geometry::disjoint(poly, other);
    }

    template <class PolyT>
    CartesianPolygon intersectionSingle(PolyT const& other) const;

    template <class PolyT>
    std::vector<CartesianPolygon> intersection(PolyT const& other) const;

    template <class PolyT>
    CartesianPolygon unionSingle(PolyT const& other) const;

    template <class PolyT>
    std::vector<CartesianPolygon> union_(PolyT const& other) const;

    template <typename MaskedImage>
    double dotProduct(MaskedImage const& image, afw::image::MaskPixel maskVal=0) const;

    BoostPolygon poly;
};

template <class PolyT>
CartesianPolygon CartesianPolygon::Impl::intersectionSingle(PolyT const& other) const
{
    std::vector<BoostPolygon> result;
    boost::geometry::intersection(poly, other, result);
    if (result.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::LengthErrorException, "Polygons have no intersection");
    }
    if (result.size() > 1) {
        throw LSST_EXCEPT(pex::exceptions::LengthErrorException,
                          (boost::format("Multiple polygons (%d) created by intersection()") %
                           result.size()).str());
    }
    return CartesianPolygon(PTR(Impl)(new Impl(result[0])));
}

template <class PolyT>
std::vector<CartesianPolygon> CartesianPolygon::Impl::intersection(PolyT const& other) const
{
    std::vector<BoostPolygon> boostResult;
    boost::geometry::intersection(poly, other, boostResult);
    std::vector<CartesianPolygon> lsstResult;
    lsstResult.reserve(boostResult.size());
    for (typename std::vector<BoostPolygon>::const_iterator i = boostResult.begin();
         i != boostResult.end(); ++i) {
        lsstResult.push_back(CartesianPolygon(PTR(Impl)(new Impl(*i))));
    }
    return lsstResult;
}

template <class PolyT>
CartesianPolygon CartesianPolygon::Impl::unionSingle(PolyT const& other) const
{
    std::vector<BoostPolygon> result;
    boost::geometry::union_(poly, other, result);
    if (result.size() != 1) {
        throw LSST_EXCEPT(pex::exceptions::LengthErrorException,
                          (boost::format("Multiple polygons (%d) created by union_()") %
                           result.size()).str());
    }
    return CartesianPolygon(PTR(Impl)(new Impl(result[0])));
}

template <class PolyT>
std::vector<CartesianPolygon> CartesianPolygon::Impl::union_(PolyT const& other) const
{
    std::vector<BoostPolygon> boostResult;
    boost::geometry::union_(poly, other, boostResult);
    std::vector<CartesianPolygon> lsstResult;
    lsstResult.reserve(boostResult.size());
    for (typename std::vector<BoostPolygon>::const_iterator i = boostResult.begin();
         i != boostResult.end(); ++i) {
        lsstResult.push_back(CartesianPolygon(PTR(Impl)(new Impl(*i))));
    }
    return lsstResult;
}

template <typename MaskedImage>
double CartesianPolygon::Impl::dotProduct(MaskedImage const& image,
                                          afw::image::MaskPixel maskVal) const
{
    afw::geom::Box2D const polyBounds = boostBoxToLsst(boost::geometry::return_envelope<BoostBox>(poly));
    afw::geom::Box2I const imageBounds = image.getBBox(afw::image::PARENT);
    int const xMin = std::max(static_cast<int>(polyBounds.getMinX()), imageBounds.getMinX());
    int const xMax = std::min(static_cast<int>(::ceil(polyBounds.getMaxX())), imageBounds.getMaxX());
    int const yMin = std::max(static_cast<int>(polyBounds.getMinY()), imageBounds.getMinY());
    int const yMax = std::min(static_cast<int>(::ceil(polyBounds.getMaxY())), imageBounds.getMaxY());
    int const x0 = image.getX0(), y0 = image.getY0();
    double sum = 0;
    for (int y = yMin; y <= yMax; ++y) {
        typename MaskedImage::Image::x_iterator i = image.getImage()->x_at(xMin - x0, y - y0);
        typename MaskedImage::Mask::x_iterator m = image.getMask()->x_at(xMin - x0, y - y0);
        for (int x = xMin; x <= xMax; ++x, ++i, ++m) {
            if (!(*m & maskVal) && boost::geometry::within(LsstPoint(x, y), poly)) {
                sum += *i;
            }
        }
    }
    return sum;
}



CartesianPolygon::CartesianPolygon(CartesianPolygon::Box const& box) :
    _impl(new CartesianPolygon::Impl(box)) {}

CartesianPolygon::CartesianPolygon(std::vector<CartesianPolygon::Point> const& vertices) :
    _impl(new CartesianPolygon::Impl(vertices)) {}

CartesianPolygon::CartesianPolygon(
    CartesianPolygon::Box const &box,
    CONST_PTR(afw::geom::XYTransform) const& transform
    ) : _impl(new CartesianPolygon::Impl())
{
    std::vector<LsstPoint> corners = boxToCorners(box);
    for (typename std::vector<LsstPoint>::iterator p = corners.begin(); p != corners.end(); ++p) {
        *p = transform->forwardTransform(*p);
    }
    boost::geometry::assign(_impl->poly, corners);
    _impl->check();
}

CartesianPolygon::CartesianPolygon(
    CartesianPolygon::Box const &box,
    afw::geom::AffineTransform const& transform
    ) : _impl(new CartesianPolygon::Impl())
{
    std::vector<LsstPoint> corners = boxToCorners(box);
    for (typename std::vector<LsstPoint>::iterator p = corners.begin(); p != corners.end(); ++p) {
        *p = transform(*p);
    }
    boost::geometry::assign(_impl->poly, corners);
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
    return _impl->overlaps(other._impl->poly);
}

bool CartesianPolygon::overlaps(Box const& box) const {
    return _impl->overlaps(box);
}

CartesianPolygon CartesianPolygon::intersectionSingle(CartesianPolygon const& other) const {
    return _impl->intersectionSingle(other._impl->poly);
}

CartesianPolygon CartesianPolygon::intersectionSingle(Box const& box) const {
    return _impl->intersectionSingle(box);
}

std::vector<CartesianPolygon> CartesianPolygon::intersection(CartesianPolygon const& other) const {
    return _impl->intersection(other._impl->poly);
}

std::vector<CartesianPolygon> CartesianPolygon::intersection(Box const& box) const {
    return _impl->intersection(box);
}

CartesianPolygon CartesianPolygon::unionSingle(CartesianPolygon const& other) const {
    return _impl->unionSingle(other._impl->poly);
}

CartesianPolygon CartesianPolygon::unionSingle(Box const& box) const {
    return _impl->unionSingle(box);
}

std::vector<CartesianPolygon> CartesianPolygon::union_(CartesianPolygon const& other) const {
    return _impl->union_(other._impl->poly);
}

std::vector<CartesianPolygon> CartesianPolygon::union_(Box const& box) const {
    return _impl->union_(box);
}

CartesianPolygon CartesianPolygon::convexHull() const
{
    BoostPolygon hull;
    boost::geometry::convex_hull(_impl->poly, hull);
    return CartesianPolygon(PTR(Impl)(new Impl(hull)));
}

PTR(afw::image::Image<float>) CartesianPolygon::createImage(afw::geom::Box2I const& bbox) const
{
    typedef afw::image::Image<float> Image;
    PTR(Image) image = boost::make_shared<Image>(bbox);
    image->setXY0(bbox.getMin());
    *image = 0.0;
    afw::geom::Box2D bounds = getBBox(); // Polygon bounds
    int xMin = std::max(static_cast<int>(bounds.getMinX()), bbox.getMinX());
    int xMax = std::min(static_cast<int>(::ceil(bounds.getMaxX())), bbox.getMaxX());
    int yMin = std::max(static_cast<int>(bounds.getMinY()), bbox.getMinY());
    int yMax = std::min(static_cast<int>(::ceil(bounds.getMaxY())), bbox.getMaxY());
    for (int y = yMin; y <= yMax; ++y) {
        double const yPixelMin = (double)y - 0.5, yPixelMax = (double)y + 0.5;
        BoostPolygon row;               // A polygon of row y
        boost::geometry::assign(row, LsstBox(afw::geom::Point2D(xMin, yPixelMin),
                                             afw::geom::Point2D(xMax, yPixelMax)));
        std::vector<BoostPolygon> intersections;
        boost::geometry::intersection(_impl->poly, row, intersections);

        if (intersections.size() == 1 && boost::geometry::num_points(intersections[0]) == 5) {
            // This row is fairly tame, and should have a long run of pixels within the polygon
            BoostPolygon const& row = intersections[0];
            std::vector<double> top, bottom;
            top.reserve(2);
            bottom.reserve(2);
            bool failed = false;
            for (std::vector<Point>::const_iterator i = row.outer().begin();
                 i != row.outer().end() - 1; ++i) {
                double const xCoord = i->getX(), yCoord = i->getY();
                if (yCoord == yPixelMin) {
                    bottom.push_back(xCoord);
                } else if (yCoord == yPixelMax) {
                    top.push_back(xCoord);
                } else {
                    failed = true;
                    break;
                }
            }
            if (!failed && top.size() == 2 && bottom.size() == 2) {
                std::sort(top.begin(), top.end());
                std::sort(bottom.begin(), bottom.end());
                int const xMin = std::min(top[0], bottom[0]);
                int const xStart = ::ceil(std::max(top[0], bottom[0])) + 1;
                int const xStop = std::min(top[1], bottom[1]) - 1;
                int const xMax = ::ceil(std::max(top[1], bottom[1]));
                pixelRowOverlap(image, _impl->poly, xMin, xStart, y);
                int x = xStart;
                for (Image::x_iterator i = image->x_at(xStart - image->getX0(), y - image->getY0());
                     x <= xStop; ++i, ++x) {
                    *i = 1.0;
                }
                pixelRowOverlap(image, _impl->poly, xStop, xMax, y);
                continue;
            }
        }

        // Last resort: do each pixel independently...
        for (typename std::vector<BoostPolygon>::const_iterator p = intersections.begin();
             p != intersections.end(); ++p) {
            double xMinRow = xMax, xMaxRow = xMin;
            std::vector<LsstPoint> const vertices = p->outer();
            for (typename std::vector<LsstPoint>::const_iterator q = vertices.begin();
                 q != vertices.end(); ++q) {
                double const x = q->getX();
                if (x < xMinRow) xMinRow = x;
                if (x > xMaxRow) xMaxRow = x;
            }

            pixelRowOverlap(image, _impl->poly, xMinRow, ::ceil(xMaxRow), y);
        }
    }
    return image;
}

double CartesianPolygon::dotProduct(
    afw::image::MaskedImage<float, afw::image::MaskPixel> const& image,
    afw::image::MaskPixel maskVal
    ) const {
    return _impl->dotProduct(image, maskVal);
}
double CartesianPolygon::dotProduct(CONST_PTR(afw::image::Image<float>) const& image) const
{
    typedef DummyMaskedImage<afw::image::Image<float>, DummyMask > MaskedImage;
    MaskedImage dummy(image, boost::make_shared<DummyMask>(), image->getBBox(afw::image::PARENT));
    return _impl->dotProduct(dummy);
}
double CartesianPolygon::dotProduct(
    CONST_PTR(afw::image::Mask<afw::image::MaskPixel>) const& mask,
    afw::image::MaskPixel maskVal
    ) const
{
    typedef afw::image::Mask<afw::image::MaskPixel> Mask;
    typedef DummyMaskedImage<DummyImage, Mask> MaskedImage;
    MaskedImage dummy(boost::make_shared<DummyImage>(), mask, mask->getBBox(afw::image::PARENT));
    return _impl->dotProduct(dummy, maskVal);
}
double CartesianPolygon::dotProduct(afw::geom::Box2I const& bbox) const
{
    typedef DummyMaskedImage<DummyImage, DummyMask> MaskedImage;
    MaskedImage dummy(boost::make_shared<DummyImage>(), boost::make_shared<DummyMask>(), bbox);
    return _impl->dotProduct(dummy);
}


}}} // namespace lsst::meas::algorithms

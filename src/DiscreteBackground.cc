#include "boost/make_shared.hpp"
#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"

#include "ndarray.h"

#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Mask.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/math/LeastSquares.h"

#include "lsst/meas/algorithms/CartesianPolygon.h"
#include "lsst/meas/algorithms/DiscreteBackground.h"

namespace {

/// Dot product of two images
///
/// Takes masked pixels into account
template <typename T1, typename T2, typename TM>
double dotProduct(CONST_PTR(lsst::afw::image::Image<T1>) const& image1,
                  CONST_PTR(lsst::afw::image::Image<T2>) const& image2,
                  CONST_PTR(lsst::afw::image::Mask<TM>) const& mask,
                  TM const maskVal             ///< Mask value of pixels to ignore
    )
{
    assert(image1->getDimensions() == image2->getDimensions());
    assert(image1->getDimensions() == mask->getDimensions());
    double sum = 0.0;

    typename lsst::afw::image::Image<T1>::const_iterator i1 = image1->begin();
    typename lsst::afw::image::Image<T2>::const_iterator i2 = image2->begin();
    typename lsst::afw::image::Mask<TM>::const_iterator m = mask->begin();
    for (; i1 != image1->end() && i2 != image2->end() && m != mask->end(); ++i1, ++i2, ++m) {
        if (!(*m & maskVal)) {
            sum += *i1 * *i2;
        }
    }
    return sum;
}

} // anonymous namespace


namespace lsst { namespace meas { namespace algorithms {

struct DiscreteBackground::Impl
{
    Impl() : polygons(), solution() {}

    Impl(Exposure const& exp, afw::geom::Box2I const& box, CONST_PTR(afw::image::Wcs) const& wcs,
         afw::image::MaskPixel const maskVal)
        : polygons(1, CartesianPolygon(afw::geom::Box2D(box),
                                       boost::make_shared<afw::image::XYTransformFromWcsPair const>(
                                           exp.getWcs(), wcs))),
          solution(),
          bbox(exp.getBBox(afw::image::PARENT)) {
        solve(exp.getMaskedImage(), maskVal);
    }
    Impl(Exposure const& exp, BoxVector const& boxList, WcsVector const& wcsList,
         afw::image::MaskPixel const maskVal)
        : polygons(), solution(), bbox(exp.getBBox(afw::image::PARENT)) {
        if (boxList.size() != wcsList.size()) {
            throw LSST_EXCEPT(pex::exceptions::LengthErrorException,
                              (boost::format("Box and Wcs length mismatch: %d vs %d") %
                               boxList.size() % wcsList.size()).str());
        }
        polygons.reserve(boxList.size());
        CONST_PTR(afw::image::Wcs) const& wcs = exp.getWcs();
        BoxVector::const_iterator b = boxList.begin();
        WcsVector::const_iterator w = wcsList.begin();
        for (int i = 0; b != boxList.end(); ++b, ++w, ++i) {
            polygons.push_back(CartesianPolygon(
                                   afw::geom::Box2D(*b),
                                   boost::make_shared<afw::image::XYTransformFromWcsPair>(wcs, *w)));
        }
        solve(exp.getMaskedImage(), maskVal);
    }
    Impl(MaskedImage const& image, PolygonVector const& _polygons, afw::image::MaskPixel maskVal)
        : polygons(), solution(), bbox(image.getBBox(afw::image::PARENT)) {
        for (PolygonVector::const_iterator i = _polygons.begin(); i != _polygons.end(); ++i) {
            if (i->overlaps(bbox)) {
                polygons.push_back(*i);
            }
        }
        solve(image, maskVal);
    }
    Impl(PolygonVector _polygons, Solution const& _solution, afw::geom::Box2I _bbox)
        : polygons(_polygons), solution(_solution), bbox(_bbox) {}

    void solve(MaskedImage const& mi, afw::image::MaskPixel const maskVal);

    void _assert() const {
        assert(polygons.size() == static_cast<size_t>(solution.getShape()[0])); // shape is signed
    }

    PolygonVector polygons;
    Solution solution;
    afw::geom::Box2I bbox;
};

//#define SOLVE_EXACT

void DiscreteBackground::Impl::solve(MaskedImage const& mi, afw::image::MaskPixel const maskVal)
{
    // Construct normal equation and solve.  Normal equation consists of matrix and vector.
    // Matrix is dot product of model and model; vector is dot product of model and data.
    size_t num = polygons.size();
    typedef ndarray::Array<double, 2, 2> Matrix;
    typedef ndarray::Array<double, 1, 1> Vector;
    Matrix matrix = ndarray::allocate(num, num);
    Vector vector = ndarray::allocate(num);

    afw::geom::Box2I const& bbox = mi.getBBox(afw::image::PARENT);
    CONST_PTR(afw::image::Image<float>) const& image = mi.getImage();
    CONST_PTR(afw::image::Mask<afw::image::MaskPixel>) const& mask = mi.getMask();

    // More worried about memory footprint than raw speed, so not caching results of Polygon.createImage().
    // XXX can do much faster dot products by making a subimage of where the polygon(s) are non-zero.
    for (size_t i = 0; i < num; ++i) {
        CartesianPolygon const iPoly = polygons[i];
#ifdef SOLVE_EXACT
        std::cerr << i << ": ";
        CONST_PTR(afw::image::Image<float>) const& iImage = iPoly.createImage(bbox);
        matrix[i][i] = dotProduct(iImage, iImage, mask, maskVal);
        vector[i] = dotProduct(image, iImage, mask, maskVal);
#else
        matrix[i][i] = iPoly.dotProduct(mask, maskVal);
        vector[i] = iPoly.dotProduct(mi, maskVal);
        if (std::isnan(vector[i])) {
            vector[i] = 0.0;
        }
#endif
        for (size_t j = i; j < num; ++j) {
            CartesianPolygon const jPoly = polygons[j];
#ifdef SOLVE_EXACT
            std::cerr << j << " ";
            CONST_PTR(afw::image::Image<float>) const& jImage = jPoly.createImage(bbox);
            matrix[i][j] = matrix[j][i] = dotProduct(iImage, jImage, mask, maskVal);
#else
            std::vector<CartesianPolygon> const intersection = iPoly.intersection(jPoly);
            double sum = 0.0;
            for (std::vector<CartesianPolygon>::const_iterator k = intersection.begin();
                 k != intersection.end(); ++k) {
                sum += k->dotProduct(mask, maskVal);
            }
            matrix[i][j] = matrix[j][i] = sum;
#endif
        }
#ifdef SOLVE_EXACT
        std::cerr << std::endl;
#endif
    }

    std::cout << matrix << std::endl;
    std::cout << vector << std::endl;

    afw::math::LeastSquares solver = afw::math::LeastSquares::fromNormalEquations(matrix, vector);
    ndarray::Array<double, 1, 1> soln1 = ndarray::copy(solver.getSolution());

    std::cout << soln1 << std::endl;
#if 1
    for (size_t i = 0; i < num; ++i) {
        if (std::isnan(soln1[i])) {
            soln1[i] = 0.0;
        }
    }
#endif
    Solution soln2 = ndarray::copy(soln1);
    solution.swap(soln2);
}

DiscreteBackground::DiscreteBackground(Exposure const& exp,
                                       afw::geom::Box2I const& box,
                                       CONST_PTR(afw::image::Wcs) const& wcs,
                                       afw::image::MaskPixel const maskVal
    ) : _impl(new Impl(exp, box, wcs, maskVal)) {}

DiscreteBackground::DiscreteBackground(Exposure const& exp,
                                       BoxVector const& boxList,
                                       WcsVector const& wcsList,
                                       afw::image::MaskPixel const maskVal
    ) : _impl(new Impl(exp, boxList, wcsList, maskVal)) {}

DiscreteBackground::DiscreteBackground(MaskedImage const& image,
                                       PolygonVector const& polygons,
                                       afw::image::MaskPixel const maskVal
    ) : _impl(new Impl(image, polygons, maskVal)) {}

DiscreteBackground::DiscreteBackground(PolygonVector const& polyList,
                                       Solution const& solution,
                                       afw::geom::Box2I const& bbox
    ) : _impl(new Impl(polyList, solution, bbox)) {}

PTR(afw::image::Image<DiscreteBackground::PixelT>) DiscreteBackground::getImage() const
{
    _impl->_assert();
    PTR(afw::image::Image<PixelT>) image;
    PolygonVector::const_iterator p = _impl->polygons.begin();
    Solution::Iterator s = _impl->solution.begin();
    for (; p != _impl->polygons.end(); ++p, ++s) {
        PTR(afw::image::Image<float>) poly = p->createImage(_impl->bbox);
        if (!image) {
            *poly *= *s;
            image = poly;
        } else {
            image->scaledPlus(*s, *poly);
        }
    }
    return image;
}

bool DiscreteBackground::operator==(DiscreteBackground const& other) const
{
    Impl const& lhs = *_impl, rhs = *other._impl;
    if (lhs.polygons.size() != rhs.polygons.size()) return false;
    lhs._assert();
    rhs._assert();
    if (lhs.bbox != rhs.bbox) return false;
    // Can't directly compare 'solution's because ndarray::Array's operator== does a shallow comparison
    if (lhs.solution.asEigen() != rhs.solution.asEigen()) return false;
    for (PolygonVector::const_iterator i = lhs.polygons.begin(), j = rhs.polygons.begin();
         i != lhs.polygons.end(); ++i, ++j) {
        if (*i != *j) return false;
    }
    return true;
}

DiscreteBackground::PolygonVector DiscreteBackground::getPolygonVector() const { return _impl->polygons; }
DiscreteBackground::Solution DiscreteBackground::getSolution() const { return _impl->solution; }
afw::geom::Box2I DiscreteBackground::getBBox() const { return _impl->bbox; }

}}} // namespace lsst::meas::algorithms

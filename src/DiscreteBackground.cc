#include "boost/make_shared.hpp"
#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/CartesianPolygon.h"
#include "lsst/meas/algorithms/DiscreteBackground.h"

namespace {

} // anonymous namespace


namespace lsst { namespace meas { namespace algorithms {

struct DiscreteBackground::Impl
{
    Impl() : polygons(), solution() {}

    Impl(Exposure const& exp, afw::geom::Box2I const& box, CONST_PTR(afw::image::Wcs) const& wcs)
        : polygons(1, CartesianPolygon(afw::geom::Box2D(box), *wcs, *exp.getWcs())), solution() { solve(exp); }
    Impl(Exposure const& exp, BoxVector const& boxList, WcsVector const& wcsList)
        : polygons(), solution() {
        if (boxList.size() != wcsList.size()) {
            throw LSST_EXCEPT(pex::exceptions::LengthErrorException,
                              (boost::format("Box and Wcs length mismatch: %d vs %d") %
                               boxList.size() % wcsList.size()).str());
        }
        polygons.reserve(boxList.size());
        afw::image::Wcs const& wcs = *exp.getWcs();
        BoxVector::const_iterator b = boxList.begin();
        WcsVector::const_iterator w = wcsList.begin();
        for (int i = 0; b != boxList.end(); ++b, ++w, ++i) {
            polygons[i] = CartesianPolygon(afw::geom::Box2D(*b), **w, wcs);
        }
        solve(exp);
    }
    Impl(PolygonVector _polygons, Solution const& _solution)
        : polygons(_polygons), solution(_solution) {}

    void solve(Exposure const& exp);

    PolygonVector polygons;
    Solution solution;
};


void DiscreteBackground::Impl::solve(Exposure const& exp)
{

}

}}} // namespace lsst::meas::algorithms

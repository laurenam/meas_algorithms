// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2008-2015 LSST/AURA.
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

#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/afw/math/Statistics.h"

#include "lsst/meas/algorithms/Variance.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/// Measurement algorithm to report the variance around an object
class VarianceAlgorithm : public Algorithm {
public:

    VarianceAlgorithm(VarianceControl const& ctrl, afw::table::Schema & schema) :
        Algorithm(ctrl),
        _resultKey(schema.addField<double>(ctrl.name, "variance at object position")),
        _maskVal(afw::image::Mask<>::getPlaneBitMask(ctrl.mask))
        {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const& exposure,
        afw::geom::Point2D const& center
        ) const;

    afw::table::Key<double> _resultKey;    // Key for accessing result
    afw::image::MaskPixel _maskVal;        // Mask bits for pixels to ignore

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(VarianceAlgorithm);
};


template <typename PixelT>
class MedianVariance : public afw::detection::FootprintFunctor<afw::image::MaskedImage<PixelT> >
{
public:
    MedianVariance(afw::image::MaskedImage<PixelT> const& image, afw::image::MaskPixel maskVal) :
        afw::detection::FootprintFunctor<afw::image::MaskedImage<PixelT> >(image),
        _maskVal(maskVal),
        _variance(),
        _stats() {
        _stats.setNanSafe(true);
    }

    virtual void reset(afw::detection::Footprint const& foot) {
        _variance.clear();
        _variance.reserve(foot.getArea());
    }

    virtual void operator()(typename afw::image::MaskedImage<PixelT>::xy_locator loc, int x, int y) {
        afw::image::MaskPixel const mask = loc.mask(0, 0);
        if (mask & _maskVal) {
            return;
        }
        _variance.push_back(loc.variance(0, 0));
    }

    PixelT getMedianVariance() const {
        return afw::math::makeStatistics(_variance, afw::math::MEDIAN, _stats).getValue();
    }

private:
    afw::image::MaskPixel _maskVal;
//    size_t _num;
    std::vector<PixelT> _variance;
    afw::math::StatisticsControl _stats;
};

template <typename PixelT>
void VarianceAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const& center
    ) const
{
    afw::geom::ellipses::Ellipse aperture(source.getShape(), center);
    aperture.scale(dynamic_cast<VarianceControl const&>(getControl()).scale);
    afw::detection::Footprint foot(aperture);
    foot.clipTo(exposure.getBBox(afw::image::PARENT));
    MedianVariance<PixelT> func(exposure.getMaskedImage(), _maskVal);
    func.apply(foot);
    source.set(_resultKey, func.getMedianVariance());
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(VarianceAlgorithm);

} // anonymous namespace


PTR(Algorithm) VarianceControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const& metadata
    ) const
{
    return boost::make_shared<VarianceAlgorithm>(*this, boost::ref(schema));
}

}}} /// namespace lsst::meas::algorithms

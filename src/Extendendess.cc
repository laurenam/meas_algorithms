// -*- lsst-c++ -*-

#include "boost/make_shared.hpp"

#include "lsst/meas/algorithms/Extendedness.h"

namespace lsst { namespace meas { namespace algorithms {

namespace {

class ExtendednessAlgorithm : public Algorithm {
public:

    ExtendednessAlgorithm(ExtendednessControl const & ctrl, afw::table::Schema & schema);

public:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

private:
    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(ExtendednessAlgorithm);

    afw::table::Key<double> _key;
};

ExtendednessAlgorithm::ExtendednessAlgorithm(
    ExtendednessControl const & ctrl, afw::table::Schema & schema
) :
    Algorithm(ctrl),
    _key(schema.addField<double>(ctrl.name, "probability of being extended"))
{}

template <typename PixelT>
void ExtendednessAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    ExtendednessControl const & ctrl = static_cast<ExtendednessControl const &>(getControl());
    source[_key] = (ctrl.fluxRatio*source.getModelFlux() + ctrl.modelErrFactor*source.getModelFluxErr())
        < (source.getPsfFlux() + ctrl.psfErrFactor*source.getPsfFluxErr()) ? 0.0 : 1.0;
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(ExtendednessAlgorithm);

} // anonymous

PTR(AlgorithmControl) ExtendednessControl::_clone() const {
    return boost::make_shared<ExtendednessControl>(*this);
}

PTR(Algorithm) ExtendednessControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<ExtendednessAlgorithm>(*this, boost::ref(schema));
}

}}} // namespace lsst::meas::algorithms

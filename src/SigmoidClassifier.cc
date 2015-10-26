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

#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/SigmoidClassifier.h"

namespace lsst { namespace meas { namespace algorithms {

namespace {

double evaluatePolynomial(
    ndarray::Array<double const,1,1> const & coefficients,
    ndarray::Vector<double,1> const & x
) {
    int const n = coefficients.template getSize<0>();
    double p = 1;
    double r = 0;
    for (int i = 0; i < n; ++i) {
        r += p*coefficients[i];
        p *= x[0];
    }
    return r;
}

template <int N>
double evaluatePolynomial(
    ndarray::Array<double const,N,N> const & coefficients,
    ndarray::Vector<double,N> const & x
) {
    double r = 0;
    int const n = coefficients.template getSize<0>();
    double p = 1.0;
    for (int i = 0; i < n; ++i) {
        r += evaluatePolynomial(coefficients[i].shallow(), x.template last<N-1>());
        p *= x[0];
    }
    return r;
}

} // anonymous


template <int N>
ndarray::Array<double const,N,N> SigmoidClassifierControl<N>::getCoefficients(
    std::string const & name
) const {
    typename CoefficientMap::const_iterator i = _coefficients.find(name);
    if (i == _coefficients.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::NotFoundException,
            (boost::format("No coefficients specified for filter '%s'") % name).str()
        );
    }
    return i->second;
}


template <int N>
void SigmoidClassifierControl<N>::setCoefficients(
    std::string const & name,
    ndarray::Array<double const,N,N> const & coefficients
) {
    _coefficients[name] = coefficients;
}


template <int N>
SigmoidClassifierAlgorithm<N>::SigmoidClassifierAlgorithm(
    SigmoidClassifierControl<N> const & ctrl,
    afw::table::Key<double> const & key
) : Algorithm(ctrl), _key(key) {}

template <int N>
template <typename PixelT>
void SigmoidClassifierAlgorithm<N>::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    ndarray::Vector<double,N> x(0.0);
    ndarray::Array<double const,N,N> coeff = getControl().getCoefficients(exposure.getFilter().getName());
    typename ndarray::Array<double const,N,N>::Index shape = coeff.getShape();
    ndarray::Vector<bool,N> active(false);
    for (int i = 0; i < N; ++i) {
        active[i] = shape[i] > 1;
    }
    getPoint(source, x, active, exposure.getCalib());
    double r = evaluatePolynomial(coeff, x);
    source.set(_key, 1.0/(1.0 + std::exp(-r)));
}

#define IMPLEMENT_ALGORITHM_INTERFACE(PIXEL) \
    template <int N> \
    void SigmoidClassifierAlgorithm<N>::_applyT( \
        lsst::afw::table::SourceRecord & source, \
        lsst::afw::image::Exposure< PIXEL > const & exposure, \
        lsst::afw::geom::Point2D const & center \
    ) const { \
        this->_apply(source, exposure, center); \
    } \
    template <int N> \
    void SigmoidClassifierAlgorithm<N>::_applyForcedT( \
        lsst::afw::table::SourceRecord & source, \
        lsst::afw::image::Exposure< PIXEL > const & exposure, \
        lsst::afw::geom::Point2D const & center, \
        lsst::afw::table::SourceRecord const & reference, \
        lsst::afw::geom::AffineTransform const & refToMeas \
    ) const { \
        this->_applyForced(source, exposure, center, reference, refToMeas); \
    }

IMPLEMENT_ALGORITHM_INTERFACE(float);
IMPLEMENT_ALGORITHM_INTERFACE(double);

#define INSTANTIATE(M) \
    template class SigmoidClassifierControl<M>; \
    template class SigmoidClassifierAlgorithm<M>

INSTANTIATE(1);
INSTANTIATE(2);
INSTANTIATE(3);

}}} // namespace lsst::meas::algorithms

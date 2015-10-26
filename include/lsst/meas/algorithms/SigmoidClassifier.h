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

#ifndef LSST_MEAS_ALGORITHMS_SigmoidClassifier_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_SigmoidClassifier_h_INCLUDED

#include <map>

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  Base class control/factory for algorithms that do star/galaxy classification
 *  by transforming a point in R^N to a scalar between 0 and 1.
 *
 *  At present, only N==1, N==2, and N==3 are supported.  In the future, this class
 *  should probably support points in arbitrary dimensions, but that would
 *  require recursive algorithms would be much more difficult to implement.
 */
template <int N>
class SigmoidClassifierControl : public AlgorithmControl {
    typedef std::map< std::string, ndarray::Array<double const,N,N> > CoefficientMap;
public:

    SigmoidClassifierControl(std::string const & name) :
        AlgorithmControl(name, 5.0)
    {}

    PTR(SigmoidClassifierControl) clone() const {
        return boost::static_pointer_cast<SigmoidClassifierControl>(_clone());
    }

    ndarray::Array<double const,N,N> getCoefficients(std::string const & name) const;

    void setCoefficients(std::string const & name, ndarray::Array<double const,N,N> const & coefficients);

private:
    CoefficientMap _coefficients;
};

template <int N>
class SigmoidClassifierAlgorithm : public Algorithm {
public:

    typedef SigmoidClassifierControl<N> Control;

    SigmoidClassifierAlgorithm(Control const & ctrl, afw::table::Key<double> const & key);

    /// Interface required by Algorithm
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

protected:

    virtual void getPoint(
        afw::table::SourceRecord & source,
        ndarray::Vector<double,N> & x,
        ndarray::Vector<bool,N> const & active,
        PTR(afw::image::Calib const) calib
    ) const = 0;

    Control const & getControl() const { return static_cast<Control const &>(getControl()); }

private:

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SigmoidClassifierAlgorithm);

    afw::table::Key<double> _key;
};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_SigmoidClassifier_h_INCLUDED

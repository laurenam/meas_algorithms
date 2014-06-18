// -*- lsst-c++ -*-
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
#ifndef LSST_MEAS_ALGORITHMS_simpleShape_h
#define LSST_MEAS_ALGORITHMS_simpleShape_h

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/algorithms/ShapeControl.h"

namespace lsst { namespace meas { namespace algorithms {

class SimpleShapeControl : public algorithms::ShapeControl {
public:

    LSST_CONTROL_FIELD(sigma, double, "Sigma of circular Gaussian used as weight function, in pixels");
    LSST_CONTROL_FIELD(nSigmaRegion, double, "Maximum radius for pixels to include, in units of sigma");

    SimpleShapeControl() :
        algorithms::ShapeControl("shape.simple"), sigma(1.5), nSigmaRegion(3) {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  Struct to hold the results of SimpleShape when we don't run it as a plugin.
 */
struct SimpleShapeResult {
    afw::geom::ellipses::Quadrupole ellipse; ///< Measured second moments.
    afw::geom::Point2D center; ///< Measured first moments, or the input center if !recentroid
    Eigen::Matrix<double,5,5> covariance; ///< Matrix of uncertainties; ordered Ixx, Iyy, Ixy, Ix, Iy.
};

class SimpleShape : public algorithms::ShapeAlgorithm {
public:

    typedef SimpleShapeControl Control;
    typedef SimpleShapeResult Result;

    SimpleShape(SimpleShapeControl const & ctrl, afw::table::Schema & schema);

    template <typename T>
    static Result measure(
        afw::geom::ellipses::Ellipse const & weight,
        afw::image::MaskedImage<T> const & image,
        double nSigmaRegion=3.0
    );

    /**
     *  @brief Convert linear raw moments into an ellipse and centroid, and return the derivative
     *         of the conversion.
     *
     *  @note This function is mainly intended for internal use, and is only exposed publically
     *        so it can be unit-tested in Python.
     *
     *  For weight function @f$w@f$ and data @f$p@f$, the "raw" moments @f$Q@f$ are defined as:
     *  @f{eqnarray*}{
     *    Q_0 &=& \sum_n w(x_n, y_n) p_n \\
     *    Q_{xx} &=& \sum_n w(x_n, y_n) x_n^2 p_n \\
     *    Q_{yy} &=& \sum_n w(x_n, y_n) y_n^2 p_n \\
     *    Q_{xy} &=& \sum_n w(x_n, y_n) x_n y_n p_n \\
     *    Q_x &=& \sum_n w(x_n, y_n) x_n p_n \\
     *    Q_y &=& \sum_n w(x_n, y_n) y_n p_n
     *  @f}
     *  whereas the converted ellipse and centroid moments are:
     *  @f{eqnarray*}{
     *    M_{xx} &=& Q_{xx} / Q_0 - Q_x^2 \\
     *    M_{xx} &=& Q_{yy} / Q_0 - Q_y^2 \\
     *    M_{xx} &=& Q_{xy} / Q_0 - Q_x Q_y \\
     *    M_x &=& Q_x / Q_0 \\
     *    M_y &=& Q_y / Q_0
     *  @f}
     *
     *  Note the slightly unusual ordering; this is for consistency with afw::geom::ellipses::Ellipse.
     */
    static Eigen::Matrix<double,5,6> convertRawMoments(
        Eigen::Matrix<double,6,1> const & q,
        afw::geom::ellipses::Quadrupole & quadrupole,
        afw::geom::Point2D & center
    );

    /**
     *  @brief Correct moments measured with a Gaussian weight function by assuming the data was also
     *         an elliptical Gaussian, and return the derivative of the correction.
     *
     *  @note This function is mainly intended for internal use, and is only exposed publically
     *        so it can be unit-tested in Python.
     *
     *  If we naively measure Gaussian-weighted moments, we'll measure the moments of the product
     *  of the weight function and the data.  What we want is the moments of the data, as if we
     *  had measured them with no weight function (but without sacrificing the S/N benefit that
     *  comes from using a weight function).  To do that, we assume the data is also an elliptical
     *  Gaussian, and "divide" the weight function from the measured moments to compute it.
     *
     *  If @f$W@f$ and @f$M@f$ are the quadruple matrices of the weight function and measurement,
     *  and @f$\eta@f$ is the measured centroid (we work in a coordinate system where the weight
     *  function is centered at the origin), then the corrected quadrupole matrix @f$C@f$ and
     *  centroid are @f$\nu@f$ are:
     *  @f{eqnarray*}{
     *    C &=& \left(M^{-1} - W^{-1}\right)^{-1} \\
     *    \nu &=& C M^{-1} \eta
     *  @f}
     */
    static Eigen::Matrix<double,5,5> correctWeightedMoments(
        afw::geom::ellipses::Quadrupole const & weight,
        afw::geom::ellipses::Quadrupole & ellipse,
        afw::geom::Point2D & center
    );

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SimpleShape);

    afw::table::KeyTuple<afw::table::Centroid> _centroidKeys;
};

}}} // namespace lsst::meas::algorithms::simpleShape

#endif // !LSST_MEAS_ALGORITHMS_simpleShape_h

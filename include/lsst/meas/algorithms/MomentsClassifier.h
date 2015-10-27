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

#ifndef LSST_MEAS_ALGORITHMS_MomentsClassifier_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_MomentsClassifier_h_INCLUDED

#include "lsst/meas/algorithms/SigmoidClassifier.h"

namespace lsst { namespace meas { namespace algorithms {


class MomentsClassifierAlgorithm;


class MomentsClassifierControl : public SigmoidClassifierControl<3> {
public:

    MomentsClassifierControl(std::string const & name);

    PTR(MomentsClassifierControl) clone() const {
        return boost::static_pointer_cast<MomentsClassifierControl>(_clone());
    }

    LSST_CONTROL_FIELD(useDeterminant, bool, "Use the determinant of the moents instead of half the trace.");

    LSST_CONTROL_FIELD(
        useSdssBeforeHsm, bool,
        "Use shape.hsm.moments only when shape.sdss fails (instead of the opposite)"
    );

    LSST_CONTROL_FIELD(flux, std::string, "Flux field to use for magnitude and/or S/N estimates");

    PTR(MomentsClassifierAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        AlgorithmMap const & others = AlgorithmMap(),
        bool isForced = false
    ) const;

protected:

    virtual PTR(AlgorithmControl) _clone() const;

    virtual PTR(Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata
    ) const;

};


class MomentsClassifierAlgorithm : public SigmoidClassifierAlgorithm<3> {
public:

    typedef MomentsClassifierControl Control;

    MomentsClassifierAlgorithm(Control const & ctrl, afw::table::Schema & schema);

protected:

    virtual void getPoint(
        afw::table::SourceRecord & source,
        ndarray::Vector<double,3> & x,
        ndarray::Vector<bool,3> const & active,
        PTR(afw::image::Calib const) calib
    ) const;

    Control const & getControl() const {
        return static_cast<Control const &>(SigmoidClassifierAlgorithm<3>::getControl());
    }

private:
    afw::table::Key< afw::table::Moments<double> > _shape1;
    afw::table::Key< afw::table::Moments<double> > _shape2;
    afw::table::Key< afw::table::Moments<double> > _shapePsf1;
    afw::table::Key< afw::table::Moments<double> > _shapePsf2;
    afw::table::Key<afw::table::Flag> _shapeFlag1;
    afw::table::Key<afw::table::Flag> _shapeFlag2;
    afw::table::Key<double> _flux;
    afw::table::Key<double> _fluxErr;
    afw::table::Key<afw::table::Flag> _fluxFlag;
    afw::table::Key<afw::table::Flag> _mainFlag;
    afw::table::Key<afw::table::Flag> _noShapeFlag;
    afw::table::Key<afw::table::Flag> _noFluxFlag;
};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_MomentsClassifier_h_INCLUDED

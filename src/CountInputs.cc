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

#include "lsst/meas/algorithms/CountInputs.h"

#include "lsst/pex/exceptions.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/// Measurement algorithm to count the number of images contributing at a point
class CountInputs : public Algorithm {
public:

    CountInputs(CountInputsControl const& ctrl, afw::table::Schema & schema) : Algorithm(ctrl) {
        _numberKey = schema.addField<int>(ctrl.name,
            "number of images contributing at center; not including any clipping");
    }

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const& exposure,
        afw::geom::Point2D const& center
        ) const {
        if (!exposure.getInfo() || !exposure.getInfo()->getCoaddInputs()) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "No coadd inputs defined");
        }
        afw::table::ExposureCatalog const& ccds = exposure.getInfo()->getCoaddInputs()->ccds;
        source.set(_numberKey, ccds.subsetContaining(center, *exposure.getWcs(), true).size());
    }

    afw::table::Key<int> _numberKey;    // Key for accessing result

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(CountInputs);
};

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(CountInputs);

} // anonymous namespace


PTR(Algorithm) CountInputsControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const& metadata
    ) const
{
    return boost::make_shared<CountInputs>(*this, boost::ref(schema));
}

}}} /// namespace lsst::meas::algorithms

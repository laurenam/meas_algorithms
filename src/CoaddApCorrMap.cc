// -*- LSST-C++ -*-

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

#include <map>
#include <vector>
#include "boost/make_shared.hpp"

#include "lsst/meas/algorithms/CoaddBoundedField.h"

#include "lsst/meas/algorithms/CoaddApCorrMap.h"

namespace lsst {
namespace meas {
namespace algorithms {

CoaddApCorrMap::CoaddApCorrMap(
    afw::table::ExposureCatalog const& catalog,
    afw::geom::Box2I const& coaddBbox,
    CONST_PTR(afw::image::Wcs) const& coaddWcs,
    std::string const & weightFieldName
    ) : afw::image::ApCorrMap()
{
    // Assemble the BoundedFields for each type
    typedef std::map<std::string const, std::vector<CoaddBoundedFieldElement> > Everything;
    Everything everything;
    afw::table::Key<double> weightKey = catalog.getSchema()[weightFieldName];
    for (afw::table::ExposureCatalog::const_iterator row = catalog.begin(); row != catalog.end(); ++row) {
        CONST_PTR(afw::image::ApCorrMap) apCorrMap = row->getApCorrMap();
        if (!apCorrMap) {
            continue;
        }
        double const weight = row->get(weightKey);
        CONST_PTR(afw::image::Wcs) wcs = row->getWcs();
        for (ApCorrMap::Iterator iter = apCorrMap->begin(); iter != apCorrMap->end(); ++iter) {
            std::string const name = iter->first;
            PTR(afw::math::BoundedField) bf = iter->second;
            everything[name].push_back(CoaddBoundedFieldElement(bf, wcs, weight));
        }
    }

    // Construct a CoaddBoundedField for each type
    for (Everything::const_iterator iter = everything.begin(); iter != everything.end(); ++iter) {
        std::string const& name = iter->first;
        std::vector<CoaddBoundedFieldElement> const& elements = iter->second;
        set(name, boost::make_shared<CoaddBoundedField>(coaddBbox, coaddWcs, elements));
    }
}

namespace {

std::string getCoaddApCorrMapPersistenceName() { return "CoaddApCorrMap"; }

} // anonymous

// Using ApCorrMap to do the persistence --- we don't add anything extra needing persistence,
// and the CoaddBoundedField has its own peristence.
CoaddApCorrMap::Factory CoaddApCorrMap::_registration(getCoaddApCorrMapPersistenceName());


std::string
CoaddApCorrMap::getPersistenceName() const
{
    return getCoaddApCorrMapPersistenceName();
}

std::string
CoaddApCorrMap::getPythonModule() const
{
    return "lsst.meas.algorithms";
}

}}} // namespace lsst::meas::algorithms



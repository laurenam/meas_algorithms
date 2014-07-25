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

#ifndef LSST_MEAS_ALGORITHMS_CoaddApCorrMap_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_CoaddApCorrMap_h_INCLUDED

#include <string>
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/ApCorrMap.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/table/Exposure.h"

namespace lsst { namespace meas { namespace algorithms {

/// Aperture corrections for a coadd
///
/// A thin veneer over ApCorrMap and CoaddBoundedField to simplify construction.
class CoaddApCorrMap : public afw::image::ApCorrMap
{
public:
    CoaddApCorrMap(
        afw::table::ExposureCatalog const& catalog,
        afw::geom::Box2I const& coaddBbox,
        CONST_PTR(afw::image::Wcs) const& coaddWcs,
        std::string const& weightFieldName = "weight"
    );

protected:
    class Factory : public afw::image::ApCorrMap::Factory {
    public:
        Factory(std::string const & name) : afw::image::ApCorrMap::Factory(name) {}
    };

private:
    static Factory _registration;

    // See afw::table::io::Persistable::getPersistenceName
    virtual std::string getPersistenceName() const;

    // See afw::table::io::Persistable::getPythonModule
    virtual std::string getPythonModule() const;

};


}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_CoaddApCorrMap_h_INCLUDED

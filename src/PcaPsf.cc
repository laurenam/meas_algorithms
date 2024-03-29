// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
/*!
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 *
 * @file
 *
 * @ingroup algorithms
 */
#include <cmath>

#include "boost/make_shared.hpp"

#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/PcaPsf.h"
#include "lsst/afw/formatters/KernelFormatter.h"
#include "lsst/afw/detection/PsfFormatter.h"
#include "lsst/meas/algorithms/KernelPsfFactory.h"

namespace afwImage = lsst::afw::image;
namespace lsst {
namespace meas {
namespace algorithms {

PcaPsf::PcaPsf(
    PTR(afw::math::LinearCombinationKernel) kernel,
    afw::geom::Point2D const & averagePosition
) : KernelPsf(kernel, averagePosition)
{
    if (!kernel) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "PcaPsf kernel must not be null");
    }
}

PTR(afw::math::LinearCombinationKernel const) PcaPsf::getKernel() const {
    return boost::static_pointer_cast<afw::math::LinearCombinationKernel const>(
        KernelPsf::getKernel()
    );
}

PTR(afw::detection::Psf) PcaPsf::clone() const {
    return boost::make_shared<PcaPsf>(*this);
}

std::string PcaPsf::getPersistenceName() const { return "PcaPsf"; }

std::string PcaPsf::getPythonModule() const { return "lsst.meas.algorithms"; }

namespace {

// registration for table persistence
KernelPsfFactory<PcaPsf,afw::math::LinearCombinationKernel> registration("PcaPsf");

} // anonymous

}}} // namespace lsst::meas::algorithms

namespace lsst { namespace afw { namespace detection {

daf::persistence::FormatterRegistration
PsfFormatter::pcaPsfRegistration = daf::persistence::FormatterRegistration(
    "PcaPsf", typeid(meas::algorithms::PcaPsf),
    lsst::afw::detection::PsfFormatter::createInstance
);

}}} // namespace lsst::afw::detection

BOOST_CLASS_EXPORT_GUID(lsst::meas::algorithms::PcaPsf, "lsst::meas::algorithms::PcaPsf")

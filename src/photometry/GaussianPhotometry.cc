// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "boost/make_shared.hpp"
#include "boost/tuple/tuple.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"
#include "lsst/meas/algorithms/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwCoord = lsst::afw::coord;

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief A class that knows how to calculate fluxes using the GAUSSIAN photometry algorithm
 * @ingroup meas/algorithms
 */
template<typename ExposureT>
class GaussianPhotometer : public Algorithm<afwDet::Photometry, ExposureT>
{
public:
    typedef Algorithm<afwDet::Photometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<GaussianPhotometer> Ptr;
    typedef boost::shared_ptr<GaussianPhotometer const> ConstPtr;

    /// Ctor
    GaussianPhotometer(double apRadius=9.0, double shiftmax=10.0, double background=0.0) :
        AlgorithmT(), _apRadius(apRadius), _shiftmax(shiftmax), _background(background) {}

    virtual std::string getName() const { return "GAUSSIAN"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<GaussianPhotometer<ExposureT> >(_apRadius, _shiftmax, _background);
    }

    virtual void configure(lsst::pex::policy::Policy const& policy) {
        if (policy.isDouble("apRadius")) {
            _apRadius = policy.getDouble("apRadius");
        } 
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        if (policy.isDouble("shiftmax")) {
            _shiftmax = policy.getDouble("shiftmax");
        }
    }

    virtual PTR(afwDet::Photometry) measureOne(ExposurePatch<ExposureT> const&, afwDet::Source const&) const;

    virtual PTR(afwDet::Photometry) measureGroups(std::vector<PTR(ExposureGroup<ExposureT>)> const&,
                                                  afwDet::Source const&) const;

private:
    double _apRadius;
    double _shiftmax;
    double _background;
};

/************************************************************************************************************/

namespace {
    
template<typename ImageT>
std::pair<double, double>
getGaussianFlux(
        ImageT const& mimage,           // the data to process
        double background,               // background level
        double xcen, double ycen,         // centre of object
        double shiftmax,                  // max allowed centroid shift
        PTR(detail::SdssShapeImpl) shape=PTR(detail::SdssShapeImpl)() // SDSS shape measurement
               )
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!shape) {
        shape = boost::make_shared<detail::SdssShapeImpl>();
    }

    if (!detail::getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, shape.get())) {
        ;                               // Should set a flag here
    } else {
        double const scale = shape->getFluxScale();
        flux = scale*shape->getI0();
        fluxErr = scale*shape->getI0Err();
    }

    return std::make_pair(flux, fluxErr);
}


/*
 * Calculate aperture correction
 *
 * The multiplier returned will correct the measured flux for an object so that if it's a PSF we'll
 * get the aperture corrected psf flux
 */
double getApertureCorrection(afwDet::Psf::ConstPtr psf, double xcen, double ycen, double shiftmax)
{
    if (!psf) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "No PSF provided for Gaussian photometry");
    }

    typedef afwDet::Psf::Image PsfImageT;
    PsfImageT::Ptr psfImage; // the image of the PSF
    PsfImageT::Ptr psfImageNoPad;   // Unpadded image of PSF
    
    int const pad = 5;
    try {
        psfImageNoPad = psf->computeImage(afwGeom::PointD(xcen, ycen));
        
        psfImage = PsfImageT::Ptr(
            new PsfImageT(psfImageNoPad->getDimensions() + afwGeom::Extent2I(2*pad))
            );
        afwGeom::BoxI middleBBox(afwGeom::Point2I(pad, pad),
                                 psfImageNoPad->getDimensions());
        
        PsfImageT::Ptr middle(new PsfImageT(*psfImage, middleBBox, afwImage::LOCAL));
        *middle <<= *psfImageNoPad;
        psfImage->setXY0(0, 0);     // SHOULD NOT BE NEEDED; psfXCen should be 0.0 and fix getGaussianFlux
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)") % xcen % ycen).str());
        throw e;
    }
    // Estimate the GaussianFlux for the Psf
    double const psfXCen = 0.5*(psfImage->getWidth() - 1); // Center of (21x21) image is (10.0, 10.0)
    double const psfYCen = 0.5*(psfImage->getHeight() - 1);
    std::pair<double, double> const result = getGaussianFlux(*psfImage, 0.0, psfXCen, psfYCen, shiftmax);
    double const psfGaussianFlux = result.first;
#if 0
    double const psfGaussianFluxErr = result.second; // NaN -- no variance in the psfImage
#endif

    // Need to correct to the PSF flux
    double psfFlux = std::accumulate(psfImageNoPad->begin(), psfImageNoPad->end(), 0.0);
    return psfFlux/psfGaussianFlux;
}
                  
} // anonymous namespace

/************************************************************************************************************/
/**
 * Calculate the desired gaussian flux
 */
template<typename ExposureT>
afwDet::Photometry::Ptr GaussianPhotometer<ExposureT>::measureOne(ExposurePatch<ExposureT> const& patch,
                                                                  afwDet::Source const& source) const
{
    typedef typename ExposureT::MaskedImageT MaskedImageT;

    CONST_PTR(ExposureT) exposure = patch.getExposure();
    CONST_PTR(afwDet::Peak) peak = patch.getPeak();

    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = peak->getFx() - mimage.getX0(); ///< object's column position in image pixel coords
    double const ycen = peak->getFy() - mimage.getY0();  ///< object's row position
    /*
     * Find the object's adaptive-moments.  N.b. it would be better to use the SdssShape measurement
     * as this code repeats the work of that measurement
     */
    std::pair<double, double> const result = getGaussianFlux(mimage, _background, xcen, ycen, _shiftmax);
    double flux = result.first;
    double fluxErr = result.second;

    /*
     * Correct the measured flux for our object so that if it's a PSF we'll
     * get the aperture corrected psf flux
     */
    double correction = getApertureCorrection(exposure->getPsf(), xcen, ycen, _shiftmax);
    flux *=    correction;
    fluxErr *= correction;

    return boost::make_shared<afwDet::Photometry>(flux, fluxErr);
}

template<typename ExposureT>
PTR(afwDet::Photometry) GaussianPhotometer<ExposureT>::measureGroups(
    std::vector<PTR(ExposureGroup<ExposureT>)> const& groups,
    afwDet::Source const& source
    ) const
{
    typedef std::vector<PTR(ExposureGroup<ExposureT>)> GroupSetT;

    PTR(afwDet::Photometry) meas = boost::make_shared<afwDet::Photometry>(); // Root node of measurements
    PTR(detail::SdssShapeImpl) masterShape = 
        boost::make_shared<detail::SdssShapeImpl>(); // Master shape for measuring fluxes; from the first image
    afwGeom::AffineTransform masterTransform; // Linear WCS for master image
    afwCoord::Coord::Ptr masterCoord; // Sky coordinates of source
    for (typename GroupSetT::const_iterator iter = groups.begin(); iter != groups.end(); ++iter) {
        CONST_PTR(ExposureGroup<ExposureT>) group = *iter;
        if (group->size() != 1) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                              "Can currently only handle one image per group");
        }

        CONST_PTR(ExposurePatch<ExposureT>) patch = (*group)[0];
        CONST_PTR(ExposureT) exposure = patch->getExposure();
        CONST_PTR(afwDet::Peak) peak = patch->getPeak();
        CONST_PTR(afwImage::Wcs) wcs = exposure->getWcs();
        typename ExposureT::MaskedImageT const& mimage = exposure->getMaskedImage();

        double const xcen = peak->getFx() - mimage.getX0(); ///< object's column position in image pixel coords
        double const ycen = peak->getFy() - mimage.getY0();  ///< object's row position

        std::pair<double, double> result;
        if (iter == groups.begin()) {
            result = getGaussianFlux(mimage, _background, xcen, ycen, _shiftmax, masterShape);
            masterCoord = wcs->pixelToSky(afwGeom::Point2D(xcen, ycen));
            masterTransform = wcs->linearizePixelToSky(masterCoord);
        } else {
            afwGeom::AffineTransform const& transform = 
                masterTransform * wcs->linearizeSkyToPixel(masterCoord);
            PTR(detail::SdssShapeImpl) shape = masterShape->transform(transform);
            result = detail::getFixedMomentsFlux(mimage, _background, xcen, ycen, *shape);
        }

        double const correction = getApertureCorrection(exposure->getPsf(), xcen, ycen, _shiftmax);
        double const flux = result.first * correction;
        double const fluxErr = result.second * correction;
        meas->add(boost::make_shared<afwDet::Photometry>(flux, fluxErr));
    }

    return meas;
}

// Declare the existence of a "GAUSSIAN" algorithm to MeasurePhotometry
LSST_DECLARE_ALGORITHM(GaussianPhotometer, afwDet::Photometry);

}}}

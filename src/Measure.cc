/// \file

#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Measure.h"

namespace image = lsst::afw::image;
namespace detection = lsst::afw::detection;
namespace algorithms = lsst::meas::algorithms;

/************************************************************************************************************/
/*
 * \brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintCentroid : public algorithms::FootprintFunctor<MaskedImageT> {
public:
    FootprintCentroid(lsst::afw::detection::Footprint const& foot, ///< The source's Footprint
                     MaskedImageT const& mimage                    ///< The image the source lives in
                     ) : algorithms::FootprintFunctor<MaskedImageT>(foot, mimage), _n(0), _sum(0), _sumx(0), _sumy(0) {}

    /// \brief Function called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

        _n++;
        _sum += val;
        _sumx += lsst::afw::image::indexToPosition(x)*val;
        _sumy += lsst::afw::image::indexToPosition(y)*val;
    }

    /// Return the number of pixels
    int getN() const { return _n; }
    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the Footprint's column centroid
    double getX() const { return this->getImage().getX0() + _sumx/_sum; }
    /// Return the Footprint's row centroid
    double getY() const { return this->getImage().getY0() + _sumy/_sum; }
private:
    int _n;
    double _sum, _sumx, _sumy;
};

/************************************************************************************************************/

template<typename MaskedImageT>
void algorithms::measureSource(lsst::afw::detection::Source::Ptr src,    ///< the Source to receive results
                               MaskedImageT& mimage,      ///< image wherein Footprint dwells
                               lsst::afw::detection::Footprint const& fp, ///< Footprint to measure
                               float background                ///< background level to subtract
                              ) {
    FootprintCentroid<MaskedImageT> centroid(fp, mimage);
    centroid.apply();
    
    src->setColc(centroid.getX());
    src->setRowc(centroid.getY());
    src->setFlux(centroid.getSum());
}

//
// Explicit instantiations
//
template void algorithms::measureSource(detection::Source::Ptr src, image::MaskedImage<float>& mimage,
                                        detection::Footprint const &fp, float background);
//
// Why do we need double images?
//
#if 1
template void algorithms::measureSource(detection::Source::Ptr src, image::MaskedImage<double>& mimage,
                                        detection::Footprint const &fp, float background);
#endif

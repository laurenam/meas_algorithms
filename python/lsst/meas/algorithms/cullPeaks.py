import numpy
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
from lsst.pipe.base import Struct

class PeakData(Struct):
    """Data about a peak

    Cached for convenience and speed.
    """
    def __init__(self, peak, maskedImage):
        center = peak.getI()
        x, y = center
        Struct.__init__(self, peak=peak, center=center, image=maskedImage.getImage().get0(x, y),
                        variance=maskedImage.getVariance().get0(x, y))


def cullPeaks(footprint, maskedImage, nSigma, polarity=True, minThreshold=0.0, display=False):
    """Culls the peaks in a footprint that are not sufficiently isolated

    For each peak, we find the highest col that you'd have to traverse to reach a still higher
    peak, and if that col's more than nSigma below your starting point, discard the peak.

    Wikipedia (https://en.wikipedia.org/wiki/Col) explains, "A col ... is a geomorphological
    term referring to the lowest point on a mountain ridge between two peaks."  It adds,
    "The height of a summit above its highest col (called the key col) is effectively a
    measure of a mountain's prominence, an important measure of the independence of its
    summit."  We are effectively measuring the prominence of our peaks, and discarding those
    that are not sufficiently prominent.

    @param footprint  The footprint for which to cull peaks
    @param maskedImage  Corresponding MaskedImage
    @param nSigma  How many sigma above local background a peak needs to be to survive
    @param polarity  Search above threshold?
    @param minThreshold  Minimum permitted col height, assumed to be below the detection threshold
    @param display  Display peaks?
    """
    if display:
        import lsst.afw.display.ds9 as ds9

    peaks = footprint.getPeaks()
    if len(peaks) <= 1:
        if display:
            ds9.dot("+", peaks[0].getIx() - maskedImage.getX0(), peaks[0].getIy() - maskedImage.getY0(),
                    ctype=ds9.GREEN, frame=1)
        return
    data = sorted((PeakData(p, maskedImage) for p in peaks), cmp=lambda x,y: cmp(x.image, y.image),
                  reverse=polarity)
    subImage = maskedImage.Factory(maskedImage, footprint.getBBox(), afwImage.PARENT, True)
    good = [data[0],]       # First peak is always good

    if display:
        ds9.dot("+", good[0].peak.getIx() - maskedImage.getX0(), good[0].peak.getIy() - maskedImage.getY0(),
                ctype=ds9.GREEN, frame=1)

    for i, peak in enumerate(data[1:], 1):
        threshold = peak.image - nSigma*numpy.sqrt(peak.variance) # Level of col
        if (not numpy.isfinite(threshold) or (polarity and threshold < minThreshold) or
            (not polarity and threshold > minThreshold)):
            if display:
                ds9.dot("x", peak.peak.getIx() - maskedImage.getX0(),
                        peak.peak.getIy() - maskedImage.getY0(), ctype=ds9.RED, frame=1)
            continue

        stops = [p.center for j, p in enumerate(data) if j != i and
                 ((polarity and p.image > threshold) or (not polarity and p.image < threshold))]
        if not afwDet.checkFootprintAtPoint(subImage.getImage(), peak.center, threshold, polarity, stops):
            good.append(peak)
            if display:
                ds9.dot("+", peak.peak.getIx() - maskedImage.getX0(), peak.peak.getIy() - maskedImage.getY0(),
                        ctype=ds9.GREEN, frame=1)
        elif display:
            ds9.dot("x", peak.peak.getIx() - maskedImage.getX0(), peak.peak.getIy() - maskedImage.getY0(),
                    ctype=ds9.YELLOW, frame=1)

    # Include only the good peaks
    peaks.clear()
    for peak in good:
        peaks.append(peak.peak)

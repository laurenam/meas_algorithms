#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests for PSF code

Run with:
   python psf.py
or
   python
   >>> import psf; psf.run()
"""

import os, sys
from math import *
import numpy
import unittest
import eups
import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as logging
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.display.ds9 as ds9
import lsst.daf.base as dafBase
import lsst.afw.display.utils as displayUtils
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.defects as defects
import lsst.meas.algorithms.utils as maUtils
import lsst.afw.cameraGeom as cameraGeom

try:
    type(verbose)
except NameError:
    verbose = 0
    logging.Trace.setVerbosity("meas.algorithms.Interp", verbose)
    logging.Trace.setVerbosity("afw.detection.Psf", verbose)
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def psfVal(ix, iy, x, y, sigma1, sigma2, b):
    """Return the value at (ix, iy) of a double Gaussian
       (N(0, sigma1^2) + b*N(0, sigma2^2))/(1 + b)
    centered at (x, y)
    """
    return (exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma1**2) +
            b*exp(-0.5*((ix - x)**2 + (iy - y)**2)/sigma2**2))/(1 + b)

class SpatialModelPsfTestCase(unittest.TestCase):
    """A test case for SpatialModelPsf"""

    @staticmethod
    def measure(footprintSet, exposure):
        """Measure a set of Footprints, returning a SourceCatalog"""
        config = measAlg.SourceMeasurementConfig()
        config.algorithms.names = ["flags.pixel", "flux.psf", "flux.naive", "flux.gaussian", "shape.sdss"]
        config.centroider.name = "centroid.sdss"
        config.algorithms["flux.naive"].radius = 3.0
        config.slots.centroid = "centroid.sdss"
        config.slots.psfFlux = "flux.psf"
        config.slots.apFlux = "flux.naive"
        config.slots.modelFlux = None
        config.slots.instFlux = None
        config.slots.shape = "shape.sdss"

        schema = afwTable.SourceTable.makeMinimalSchema()
        measureSources = config.makeMeasureSources(schema)
        catalog = afwTable.SourceCatalog(schema)
        config.slots.setupTable(catalog.table)

        if display:
            ds9.mtv(exposure)

        footprintSet.makeSources(catalog)

        for i, source in enumerate(catalog):
            measureSources.applyWithPeak(source, exposure)

        return catalog

    def setUp(self):
        width, height = 110, 301

        self.mi = afwImage.MaskedImageF(afwGeom.ExtentI(width, height))
        self.mi.set(0)
        sd = 3                          # standard deviation of image
        self.mi.getVariance().set(sd*sd)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.FWHM = 5
        self.ksize = 31                      # size of desired kernel

        sigma1 = 1.75
        sigma2 = 2*sigma1

        self.exposure = afwImage.makeExposure(self.mi)
        self.exposure.setPsf(measAlg.DoubleGaussianPsf(self.ksize, self.ksize,
                                                    1.5*sigma1, 1, 0.1))
        self.exposure.setDetector(cameraGeom.Detector(cameraGeom.Id(1), False, 1.0))
        self.exposure.getDetector().setDistortion(None) #cameraGeom.NullDistortion())
        
        #
        # Make a kernel with the exactly correct basis functions.  Useful for debugging
        #
        basisKernelList = afwMath.KernelList()
        for sigma in (sigma1, sigma2):
            basisKernel = afwMath.AnalyticKernel(self.ksize, self.ksize,
                                                 afwMath.GaussianFunction2D(sigma, sigma))
            basisImage = afwImage.ImageD(basisKernel.getDimensions())
            basisKernel.computeImage(basisImage, True)
            basisImage /= numpy.sum(basisImage.getArray())

            if sigma == sigma1:
                basisImage0 = basisImage
            else:
                basisImage -= basisImage0

            basisKernelList.append(afwMath.FixedKernel(basisImage))

        order = 1                                # 1 => up to linear
        spFunc = afwMath.PolynomialFunction2D(order)

        exactKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
        exactKernel.setSpatialParameters([[1.0, 0,          0],
                                          [0.0, 0.5*1e-2, 0.2e-2]])
        self.exactPsf = measAlg.PcaPsf(exactKernel)        

        rand = afwMath.Random()               # make these tests repeatable by setting seed

        addNoise = True

        if addNoise:
            im = self.mi.getImage()
            afwMath.randomGaussianImage(im, rand) # N(0, 1)
            im *= sd                              # N(0, sd^2)
            del im

        xarr, yarr = [], []

        for x, y in [(20, 20), (60, 20), 
                     (30, 35),
                     (50, 50),
                     (20, 90), (70, 160), (25, 265), (75, 275), (85, 30),
                     (50, 120), (70, 80),
                     (60, 210), (20, 210),
                     ]:
            xarr.append(x)
            yarr.append(y)

        for x, y in zip(xarr, yarr):
            dx = rand.uniform() - 0.5   # random (centered) offsets
            dy = rand.uniform() - 0.5

            k = exactKernel.getSpatialFunction(1)(x, y) # functional variation of Kernel ...
            b = (k*sigma1**2/((1 - k)*sigma2**2))       # ... converted double Gaussian's "b"

            #flux = 80000 - 20*x - 10*(y/float(height))**2
            flux = 80000*(1 + 0.1*(rand.uniform() - 0.5))
            I0 = flux*(1 + b)/(2*numpy.pi*(sigma1**2 + b*sigma2**2))
            for iy in range(y - self.ksize//2, y + self.ksize//2 + 1):
                if iy < 0 or iy >= self.mi.getHeight():
                    continue

                for ix in range(x - self.ksize//2, x + self.ksize//2 + 1):
                    if ix < 0 or ix >= self.mi.getWidth():
                        continue

                    I = I0*psfVal(ix, iy, x + dx, y + dy, sigma1, sigma2, b)
                    Isample = rand.poisson(I) if addNoise else I
                    self.mi.getImage().set(ix, iy, self.mi.getImage().get(ix, iy) + Isample)
                    self.mi.getVariance().set(ix, iy, self.mi.getVariance().get(ix, iy) + I)
        # 
        bbox = afwGeom.BoxI(afwGeom.PointI(0,0), afwGeom.ExtentI(width, height))
        self.cellSet = afwMath.SpatialCellSet(bbox, 100)

        self.footprintSet = afwDetection.FootprintSet(self.mi, afwDetection.Threshold(100), "DETECTED")

        self.catalog = SpatialModelPsfTestCase.measure(self.footprintSet, self.exposure)

        for source in self.catalog:
            try:
                cand = measAlg.makePsfCandidate(source, self.exposure)
                self.cellSet.insertCandidate(cand)

            except Exception, e:
                print e
                continue

    def tearDown(self):
        del self.cellSet
        del self.exposure
        del self.mi
        del self.exactPsf
        del self.footprintSet
        del self.catalog

    @staticmethod
    def setupDeterminer(exposure, nEigenComponents=3, starSelectorAlg="secondMoment"):
        """Setup the starSelector and psfDeterminer"""
        starSelectorFactory = measAlg.starSelectorRegistry[starSelectorAlg]
        starSelectorConfig = starSelectorFactory.ConfigClass()

        if starSelectorAlg == "secondMoment":
            starSelectorConfig.clumpNSigma = 5.0
            starSelectorConfig.histSize = 14
            starSelectorConfig.badFlags = ["flags.pixel.edge",
                                           "flags.pixel.interpolated.center",
                                           "flags.pixel.saturated.center",
                                           "flags.pixel.cr.center",
                                           ]
        elif starSelectorAlg == "objectSize":
            starSelectorConfig.sourceFluxField = "flux.gaussian"
            starSelectorConfig.badFlags = ["flags.pixel.edge",
                                           "flags.pixel.interpolated.center",
                                           "flags.pixel.saturated.center",
                                           "flags.pixel.cr.center",
                                           ]
            starSelectorConfig.widthStdAllowed = 0.5
            
        starSelector = starSelectorFactory(starSelectorConfig)
        
        psfDeterminerFactory = measAlg.psfDeterminerRegistry["pca"]
        psfDeterminerConfig = psfDeterminerFactory.ConfigClass()
        width, height = exposure.getMaskedImage().getDimensions()
        psfDeterminerConfig.sizeCellX = width
        psfDeterminerConfig.sizeCellY = height//3
        psfDeterminerConfig.nEigenComponents = nEigenComponents
        psfDeterminerConfig.spatialOrder = 1
        psfDeterminerConfig.kernelSizeMin = 31
        psfDeterminerConfig.nStarPerCell = 0
        psfDeterminerConfig.nStarPerCellSpatialFit = 0 # unlimited
        psfDeterminer = psfDeterminerFactory(psfDeterminerConfig)

        return starSelector, psfDeterminer


    def subtractStars(self, exposure, catalog, chi_lim=-1):
        """Subtract the exposure's PSF from all the sources in catalog"""
        mi, psf = exposure.getMaskedImage(), exposure.getPsf()

        subtracted =  mi.Factory(mi, True)

        for s in catalog:
            xc, yc = s.getX(), s.getY()
            bbox = subtracted.getBBox(afwImage.PARENT)
            if bbox.contains(afwGeom.PointI(int(xc), int(yc))):
                try:
                    measAlg.subtractPsf(psf, subtracted, xc, yc)
                except:
                    pass

        chi = subtracted.Factory(subtracted, True)
        var = subtracted.getVariance()
        numpy.sqrt(var.getArray(), var.getArray()) # inplace sqrt
        chi /= var

        if display:
            ds9.mtv(subtracted, title="Subtracted", frame=1)
            ds9.mtv(chi, title="Chi", frame=2)
            ds9.mtv(psf.computeImage(afwGeom.Point2D(xc, yc)), title="Psf", frame=3)
            ds9.mtv(mi, frame=4, title="orig")
            kern = psf.getKernel()
            kimg = afwImage.ImageD(kern.getWidth(), kern.getHeight())
            kern.computeImage(kimg, True, xc, yc)
            ds9.mtv(kimg, title="kernel", frame=5)

        chi_min, chi_max = numpy.min(chi.getImage().getArray()),  numpy.max(chi.getImage().getArray())
        if False:
            print chi_min, chi_max

        if chi_lim > 0:
            self.assertGreater(chi_min, -chi_lim)
            self.assertLess(   chi_max,  chi_lim)

    def testPsfDeterminer(self):
        """Test the (PCA) psfDeterminer"""

        for starSelectorAlg in ["secondMoment",
                                "objectSize",
                                ]:
            print "Using %s star selector" % (starSelectorAlg)

            starSelector, psfDeterminer = \
                SpatialModelPsfTestCase.setupDeterminer(self.exposure,
                                                        nEigenComponents=2, starSelectorAlg=starSelectorAlg)
            metadata = dafBase.PropertyList()
            psfCandidateList = starSelector.selectStars(self.exposure, self.catalog)
            psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)
            self.exposure.setPsf(psf)

            chi_lim = 5.0
            self.subtractStars(self.exposure, self.catalog, chi_lim)

    def testPsfDeterminerSubimage(self):
        """Test the (PCA) psfDeterminer on subImages"""

        w, h = self.exposure.getDimensions()
        x0, y0 = int(0.35*w), int(0.45*h)
        bbox = afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(w - x0, h - y0))
        subExp = self.exposure.Factory(self.exposure, bbox)

        starSelector, psfDeterminer = \
            SpatialModelPsfTestCase.setupDeterminer(subExp, nEigenComponents=2,
                                                    starSelectorAlg="objectSize")
        metadata = dafBase.PropertyList()
        #
        # Only keep the sources that lie within the subregion (avoiding lots of log messages)
        #
        def trimCatalogToImage(exp, catalog):
            trimmedCatalog = afwTable.SourceCatalog(catalog.table.clone())
            for s in catalog:
                if exp.getBBox(afwImage.PARENT).contains(afwGeom.PointI(s.getCentroid())):
                    trimmedCatalog.append(trimmedCatalog.table.copyRecord(s))

            return trimmedCatalog

        psfCandidateList = starSelector.selectStars(subExp,
                                                    trimCatalogToImage(subExp, self.catalog))
        psf, cellSet = psfDeterminer.determinePsf(subExp, psfCandidateList, metadata)
        subExp.setPsf(psf)

        # Test how well we can subtract the PSF model.  N.b. using self.exposure is an extrapolation
        for exp, chi_lim in [(subExp, 4.5),
                             (self.exposure.Factory(self.exposure,
                                                    afwGeom.BoxI(afwGeom.PointI(0, 100),
                                                                 (afwGeom.PointI(w-1, h-1)))), 7.5),
                             (self.exposure, 19),
                             ]:
            cat = trimCatalogToImage(exp, self.catalog)
            exp.setPsf(psf)
            self.subtractStars(exp, cat, chi_lim)

    def testPsfDeterminerNEigen(self):
        """Test the (PCA) psfDeterminer when you ask for more components than acceptable stars"""

        starSelector, psfDeterminer = SpatialModelPsfTestCase.setupDeterminer(self.exposure,
                                                                              nEigenComponents=3)
        metadata = dafBase.PropertyList()
        psfCandidateList = starSelector.selectStars(self.exposure, self.catalog)
        psfCandidateList, nEigen = psfCandidateList[0:4], 2 # only enough stars for 2 eigen-components
        psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)

        self.assertEqual(psf.getKernel().getNKernelParameters(), nEigen)

    def testCandidateList(self):
        self.assertFalse(self.cellSet.getCellList()[0].empty())
        self.assertTrue(self.cellSet.getCellList()[1].empty())
        self.assertFalse(self.cellSet.getCellList()[2].empty())
        self.assertTrue(self.cellSet.getCellList()[3].empty())

        stamps = []
        stampInfo = []
        for cell in self.cellSet.getCellList():
            for cand in cell:
                #
                # Swig doesn't know that we inherited from SpatialCellMaskedImageCandidate;  all
                # it knows is that we have a SpatialCellCandidate, and SpatialCellCandidates
                # don't know about getMaskedImage;  so cast the pointer to SpatialCellMaskedImageCandidate<float>
                # and all will be well
                #
                cand = afwMath.cast_SpatialCellMaskedImageCandidateF(cell[0])
                width, height = 29, 25
                cand.setWidth(width); cand.setHeight(height);

                im = cand.getMaskedImage()
                stamps.append(im)

                self.assertEqual(im.getWidth(), width)
                self.assertEqual(im.getHeight(), height)
        
        if False and display:
            mos = displayUtils.Mosaic()
            mos.makeMosaic(stamps, frame=2)

    def testRejectBlends(self):
        """Test the PcaPsfDeterminer blend removal

        We give it a single blended source, asking it to remove blends,
        and check that it barfs in the expected way.
        """

        factory = measAlg.psfDeterminerRegistry["pca"]
        config = factory.ConfigClass()
        config.doRejectBlends = True
        psfDeterminer = factory(config)

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("position", afwGeom.Point2D, doc="Position")
        catalog = afwTable.SourceCatalog(schema)
        catalog.defineCentroid("position")
        source = catalog.addNew()

        foot = afwDetection.Footprint(afwGeom.Point2I(123, 45), 6, self.exposure.getBBox())
        peaks = foot.getPeaks()
        peaks.push_back(afwDetection.Peak(123, 45, 6))
        peaks.push_back(afwDetection.Peak(126, 47, 5))
        source.setFootprint(foot)

        candidates = [measAlg.makePsfCandidate(source, self.exposure)]
        metadata = dafBase.PropertyList()

        with self.assertRaises(RuntimeError) as cm:
            psfDeterminer.determinePsf(self.exposure, candidates, metadata)
        self.assertEqual(cm.exception.message, "All PSF candidates removed as blends")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SpatialModelPsfTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the utilsTests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

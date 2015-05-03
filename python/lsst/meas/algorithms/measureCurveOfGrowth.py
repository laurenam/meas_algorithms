#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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
import sys
import numpy as np

import lsstDebug
import lsst.pex.config as pex_config
import lsst.pipe.base as pipe_base
import lsst.afw.display.ds9 as ds9
from .algorithmsLib import getApCorrRegistry

try:
    import matplotlib.pyplot as plt
    fig = None
except ImportError as e:
    print e
    plt = None

__all__ = ("CurveOfGrowthMeasurementConfig", "CurveOfGrowthMeasurementTask", "CurveOfGrowth")

class CurveOfGrowthMeasurementConfig(pex_config.Config):
    nAperture = pex_config.Field(
        doc = "Maximum number of aperture fluxes to use (all, if None)",
        dtype = int,
        default = None
        )
    badFlags = pex_config.ListField(
        doc = """List of flags which cause a source to be rejected as bad

        N.b. may contain globs such as *.flags.pixel.edge""",
        dtype = str,
        default = ["*.flags.pixel.edge",
               ]
        )
    classificationMax = pex_config.Field(
        doc = "Maximum value of classification.extendedness for an object to be included in curve of growth",
        dtype = float,
        default = np.nan,
    )
    psfFluxMin = pex_config.Field(
        doc = "Minimum value of psf.flux for an object to be included in curve of growth",
        dtype = float,
        default = 5e4,
        check = lambda x: x >= 0.0,
    )
    nSaturated =  pex_config.Field(
        doc = "How many of the brightest saturated candidates to use (0 => all)",
        dtype = int,
        default = 500,
        check = lambda x: x >= 0,
    )
    nNonSaturated =  pex_config.Field(
        doc = "How many of the brightest non-saturated candidates to use (0 => all)",
        dtype = int,
        default = 200,
        check = lambda x: x >= 0,
    )
    skyNoiseFloor =  pex_config.Field(
        doc = "The per-pixel std. dev. of our knowledge of the sky, to be added in quadrature",
        dtype = float,
        default = 0.0,
        check = lambda x: x >= 0,
    )
    maxRChi2 = pex_config.ListField(
        doc = """List of values of reduced chi^2 that should be applied in order to clip sources

        N.b. these values are so large because of contamination in the annuli by to e.g. the faint wings of
        neighbouring objects.  The real solution here is to write cleverer aperture flux code, \'a la SDSS
        """,
        dtype = float,
        default = [10000, 1000],
        itemCheck = lambda x: x > 0,
    )
    finalEstimationAlgorithm = pex_config.ChoiceField(
        doc="""After perfoming an ML estimation of the curve of growth to estimate the relative
        weights we may undertake a per-radial-point re-estimation of the curve of growth;
        finalEstimationAlgorithm dictates how to estimate the value of each point.
        """,
        dtype=str,
        default = "mean",
        allowed={
            "mean": "weighted mean",
            "median": "median",
            None: "no clipping",
            }
        )
    fracInterpolatedMax = pex_config.Field(
        doc = """Maximum fraction of the pixels in an annulus which may be labelled INTRP

        Annuli with more than fracInterpolatedMax interpolated pixels are rejected
        """,
        dtype = float,
        default = 0.1,
    )
    minAnnularFlux = pex_config.Field(
        doc = """Minimum flux measured in an annulus to be included in the curve of growth

        Annuli with fluxes of less than minAnnularFlux are rejected
        """,
        dtype = float,
        default = 0,
    )
    def validate(self):
        pex_config.Config.validate(self)
        # Allow None to be used as an equivalent for "None"
        if self.finalEstimationAlgorithm is "NONE":
            self.finalEstimationAlgorithm = None

## \addtogroup LSST_task_documentation
## \{
## \page curveOfGrowthMeasurementTask
## \ref CurveOfGrowthMeasurementTask_ "CurveOfGrowthMeasurementTask"
## \copybrief CurveOfGrowthMeasurementTask
## \}

class CurveOfGrowthMeasurementTask(pipe_base.Task):
    """!
\anchor CurveOfGrowthMeasurementTask_
\brief Estimate a curve of growth from the aperture fluxes of a set of Sources

The measured values of aperture.flux define points on the cumulative radial profile for an object
over a range of radii, but
    - The innermost values may be missing if the star is saturated
    - The outermost values may have too low a signal-to-noise to be useful
It is often very useful to have measurements that spans the full available radial range, and
this can be created by combining the aperture fluxes of a number of objects.

\section meas_algorithms_curveOfGrowthMeasurement_Contents Contents

 - \ref meas_algorithms_curveOfGrowthMeasurement_Purpose
 - \ref meas_algorithms_curveOfGrowthMeasurement_Initialize
 - \ref meas_algorithms_curveOfGrowthMeasurement_Invoke
 - \ref meas_algorithms_curveOfGrowthMeasurement_Config
 - \ref meas_algorithms_curveOfGrowthMeasurement_Debug
 - \ref meas_algorithms_curveOfGrowthMeasurement_Example

\section meas_algorithms_curveOfGrowthMeasurement_Purpose    Description

\copybrief CurveOfGrowthMeasurementTask

This task builds a curve of growth using overlapping portions of the radial profiles of multiple
sources.  No one source is required to have measurements for all parts of the profile;  in particular
the code is intended to piece together faint objects (with well measured centres and poorly measured
wings) and very bright objects (for which the center is saturated, but the wing well measured).

The basic algorithm is a maximum likelihood fit, described <a href="curveOfGrowth.pdf">here</a>.
For HSC (and probably LSST) it proved necessary to iterate, clipping objects that don't fit an
initial robust estimate of the profile.  See the configuration parameter
\link CurveOfGrowthMeasurementConfig.maxRChi2\endlink, which gives the list of thresholds in reduced chi^2
used in the clipping.

\section meas_algorithms_curveOfGrowthMeasurement_Initialize Task initialisation

\copydoc init

\section meas_algorithms_curveOfGrowthMeasurement_Invoke     Invoking the Task

\copydoc run

\section meas_algorithms_curveOfGrowthMeasurement_Config       Configuration parameters

See \ref CurveOfGrowthMeasurementConfig

\section meas_algorithms_curveOfGrowthMeasurement_Debug              Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c --debug to import \b debug.py from your \c PYTHONPATH; see \ref baseDebug for more about \b debug.py files.

The available variables in CurveOfGrowthMeasurementTask are:
- \c display If False, disable all plotting.  If true, 
 - if displayImage is non-zero, display the exposure on ds9's on frame int(displayImage).
     Annuli selected as inputs to the curve of growth are shown (yellow for good points, blue for bad).
  - If displayImageIds is true, sources rejected are labelled by their objectIds drawn in red,
    those accepted are drawn in green
 - If plotProfiles is True, use matplotlib to show the input and final curves of growth
      using \link CurveOfGrowth.plot\endlink.

\em N.b. the exposure isn't actually displayed (the task doesn't have access to the exposure).  See
    \link meas_algorithms_curveOfGrowthMeasurement_Example\endlink for a workaround.

\section meas_algorithms_curveOfGrowthMeasurement_Example    A complete example of using CurveOfGrowthMeasurementTask

This code is in \link measAlgTasks.py\endlink in the examples directory, and can be run as \em e.g.
\code
examples/measAlgTasks.py --ds9
\endcode

See the description in \link meas_algorithms_measurement_Example\endlink

<HR>
To investigate the \ref meas_algorithms_curveOfGrowthMeasurement_Debug, put something like
\code{.py}
import lsstDebug
def DebugInfo(name):
    di = lsstDebug.getInfo(name)        # N.b. lsstDebug.Info(name) would call us recursively
    displayImage = True
    if name == "lsst.meas.algorithms.measureCurveOfGrowth":
        di.display = True
        di.displayImage = displayImage
        di.displayImageIds = True
        di.plotProfiles= True
        di.normalize = True
        di.showRadialProfile = False

    if name == "lsst.meas.algorithms.detection":
        di.display = displayImage

    return di

lsstDebug.Info = DebugInfo
\endcode
into your debug.py file and run measAlgTasks.py with the \c --debug flag.

The use of \c displayImage is a workaround for the inability of CurveOfGrowthMeasurementTask to display
the exposure.

    """

    ConfigClass = CurveOfGrowthMeasurementConfig
    _DefaultName = "curveOfGrowthMeasurement"

    def init(self, schema=None, **kwds):
        """! Create the task to estimate curves of growth

        \param schema If non-None, the measurement schema.  It will be augmented by entries:
        <DL>
        <DT> curveOfGrowth.candidate
        <DD> set if source was considered for use in constructing the curve of growth
        <DT> curveOfGrowth.used
        <DD> set if source was used in constructing the curve of growth
        (n.b. you are responsible for ensuring that they are added to the output table,
        usually by creating the table using the schema that was passed to this task)
        </DL>

        \param kwds Other arguments (including the Config) which are passed to the base class
        """
        self.__init__(schema, **kwds)

    def __init__(self, schema=None, **kwds):
        pipe_base.Task.__init__(self, **kwds)

        if schema:
            self.curveOfGrowthCandidateKey = schema.addField(
                "curveOfGrowth.candidate", type="Flag",
                doc="set if source was considered for use in constructing the curve of growth"
            )
            self.curveOfGrowthUsedKey = schema.addField(
                "curveOfGrowth.used", type="Flag",
                doc="set if source was used in constructing the curve of growth"
            )
        else:
            self.curveOfGrowthCandidateKey = None
            self.curveOfGrowthUsedKey = None

    def run(self, *catalogList):
        """!Estimate the curve of growth from a set of Sources by invoking run(catalog)
        \param catalogList  A list of SourceCatalogs

        The Sources must have these fields set:
         - flux.aperture (and other fields set by the flux.aperture algorithm)
         - all fields listed in \link CurveOfGrowthMeasurementConfig.badFlags\endlink
         - classification.extendedness
        and may have deblend.nchild set.


        \return a Struct with:
           - curveOfGrowth a CurveOfGrowth object
        """
        #
        # Decide which objects to use.  Because the catalog may not be contiguous
        # we can't numpy views into columns, so define the needed keys and perform
        # the checks on each source
        #
        sch = None
        for catalog in catalogList:
            if sch is None:
                sch = catalog.getSchema()
            else:
                assert sch == catalog.getSchema(), "Schema mismatch"

        try:
            deblend_nchild = sch.find("deblend.nchild").key
        except LookupError:
            deblend_nchild = None
        # badFlags may contain globs such as *.flags.pixel.edge which extract expands
        badFlags = []
        for badFlag in self.config.badFlags:
            badFlags += [sch.find(f).key for f in sch.extract(badFlag)]
                
        classification_extendedness = sch.find("classification.extendedness").key
        #
        # Now we've chosen our sources, add them to the nascent curve of growth
        cog = CurveOfGrowth(self.curveOfGrowthCandidateKey,
                            self.curveOfGrowthUsedKey,
                            self.config.fracInterpolatedMax, self.config.minAnnularFlux,
                            nAperture=self.config.nAperture,
                            skyNoiseFloor=self.config.skyNoiseFloor
        )

        for catalog in catalogList:
            for s in catalog:
                if (deblend_nchild is not None and s.get(deblend_nchild) > 0) or \
                   s.getPsfFlux() < self.config.psfFluxMin or \
                   (np.isfinite(self.config.classificationMax) and
                    s.get(classification_extendedness) > self.config.classificationMax):
                    continue

                for badFlag in badFlags:
                    if s.get(badFlag):
                        continue

                cog.addSource(s)
        #
        # and actually estimate it
        #
        cog.estimate(maxRChi2=self.config.maxRChi2,
                     nSaturated=self.config.nSaturated,
                     nNonSaturated=self.config.nNonSaturated,
                     finalEstimationAlgorithm=self.config.finalEstimationAlgorithm)
        #
        if lsstDebug.Info(__name__).display:
            if lsstDebug.Info(__name__).displayImage:
                
                frame = lsstDebug.Info(__name__).displayImage
                frame = 0 if frame is True else int(frame)
                if False:               # we don't actually have an exposure to display
                    ds9.mtv(exposure, frame=frame)
                else:
                    ds9.erase(frame=frame)

                for catalog in catalogList:
                    if catalog.getMetadata().exists("flux_aperture_radii"):
                        radii = catalog.getMetadata().get("flux_aperture_radii")

                        with ds9.Buffering():
                            for s in catalog:
                                if not s.get("curveOfGrowth.candidate"):
                                    continue

                                xy = s.getCentroid()
                                nInterpPixel = s.get("flux.aperture.nInterpolatedPixel")
                                nInterpOld = 0    # number of interpolated pixels inside previous radius
                                radiusOld = 0

                                for i in range(s.get("flux.aperture.nProfile")):
                                    nInterpAperture = nInterpPixel[i] - nInterpOld
                                    area = np.pi*(radii[i]**2 - radiusOld**2)
                                    nInterpOld, radiusOld = nInterpPixel[i], radii[i]

                                    ctype=ds9.YELLOW if \
                                           nInterpAperture/area < self.config.fracInterpolatedMax else ds9.BLUE

                                    ds9.dot('o', *xy, size=radii[i], frame=frame, ctype=ctype)

                                if lsstDebug.Info(__name__).displayImageIds:
                                    ds9.dot("%d" % (s.getId()),
                                            xy[0] + 0.85*radii[i], xy[1] + 0.85*radii[i], frame=frame,
                                            ctype=ds9.GREEN if s.get("curveOfGrowth.used") else ds9.RED)

            if lsstDebug.Info(__name__).plotProfiles:
                fig = cog.plot(normalize=lsstDebug.Info(__name__).normalize,
                                 showRadialProfile=lsstDebug.Info(__name__).showRadialProfile)
                fig.show()

                while True:
                    try:
                        reply = raw_input("continue? [c h(elp) q(uit) p(db)] ").strip()
                    except EOFError:
                        reply = None
                    if not reply:
                        reply = "c"

                    if reply:
                        if reply[0] == "h":
                            print """\
        At this prompt, you can continue with almost any key; 'p' enters pdb, and 'h' prints this text
        """
                        elif reply[0] == "p":
                            import pdb; pdb.set_trace()
                        elif reply[0] == 'q':
                            sys.exit(1)
                        else:
                            break

        return pipe_base.Struct(
            curveOfGrowth = cog
            )

class CurveOfGrowthResult(object):
    """!Result of the Curve of Growth calculation

    The CurveOfGrowth class manages the calculation, while this class manages
    the result and its application. The main reason for the separation is that
    the CurveOfGrowthResult is picklable, so it can be stored or transferred.
    """

    MATCHES_FIRST = 1 # apply() only to first element of matches
    MATCHES_SECOND = 2 # apply() only to second element of matches
    MATCHES_BOTH = 3 # apply() to both elements of matches

    def __init__(self, apertureFlux, apertureFluxErr):
        """!Constructor

        \param apertureFlux  ndarray of aperture flux values
        \param apertureFluxErr  ndarray of aperture flux error values
        """
        self.apertureFlux = apertureFlux
        self.apertureFluxErr = apertureFluxErr

    def getRatio(self, inner, outer):
        """!Return the flux ratio between two radii and the error on the ratio.

        \param inner     Index of the inner radius (numerator of the ratio); if -ve relative to end of array
        \param outer     Index of the outer radius (denominator of the ratio); if -ve relative to end of array
        """
        if inner < 0:
            inner += len(self.apertureFlux)
        if outer < 0:
            outer += len(self.apertureFlux)

        if outer < inner:
            raise ValueError("Inner index (%s) is larger than outer index (%d)" % (inner, outer))
        fInner = self.apertureFlux[inner]
        fOuter = self.apertureFlux[outer]
        ratio = fInner / fOuter
        # To compute error on the ratio, we propagate errors on:
        #  ratio = fInner / (fInner + delta)
        # where delta = fOuter - fInner
        # because the errors on f_inner and delta are independent,
        # while the errors on f_outer and f_inner are not.
        fInnerVar = self.apertureFluxErr[inner]**2
        delta = fOuter - fInner
        deltaVar = self.apertureFluxErr[outer]**2 - fInnerVar
        ratioErr = (fInnerVar*delta**2 + deltaVar*fInner**2)**0.5 / fOuter**2
        return ratio, ratioErr

    def apply(self, measurementConfig, catalog=None, matches=None, algorithms=None, matchesType=None,
              calib=None, apCorr=None, log=None):
        """!Apply the results of the curve of growth analysis to a source catalog

        We deduce the correction that needs to be applied from inspecting the
        'measurementConfig' (a SourceMeasurementConfig) and finding the aperture
        being used for calibration.  The correction is then applied to the flux
        of each of the specified 'algorithms'. If 'algorithms' is None, the set
        of algorithms in the aperture correction registry is used.

        Note that we do *not* correct the errors in the fluxes, since they are a
        systematic error in the zero point, rather than errors in the individual
        fluxes.

        If a 'catalog' is provided, the following entries will be added to its
        metadata:
          * CURVE_OF_GROWTH.CORRECTION.VALUE: the correction
          * CURVE_OF_GROWTH.CORRECTION.ERROR: the error in the correction
          * CURVE_OF_GROWTH.CORRECTION.RADIUS.FROM: inner radius
          * CURVE_OF_GROWTH.CORRECTION.RADIUS.TO: outer radius
          * CURVE_OF_GROWTH.CORRECTION.ALGORITHMS: list of algorithms corrected

        If you provide 'matches' to be corrected, you may also specify which part
        of the matches will be corrected through the 'matchesType' parameter:
          * CurveOfGrowthResult.MATCHES_FIRST: correct only first elements
          * CurveOfGrowthResult.MATCHES_SECOND: correct only second elements
          * CurveOfGrowthResult.MATCHES_BOTH: correct both first and second elements
        If this parameter is None, it defaults to MATCHES_BOTH.

        \param measurementConfig    Configuration for source measurement
        \param catalog    Source catalog to have flux measurements corrected
        \param matches    Source matches to be corrected
        \param algorithms    Iterable of algorithms to correct, or None
        \param matchesType    Which element of matches to correct, or None for both
        \param calib    Calib object to be corrected
        \param apCorr    Aperture corrections to be corrected
        \param log    Log object for logging the factor that's applied
        """
        # Figure out which aperture corresponds to our calibration aperture
        # This requires assuming a parameter name for the aperture;
        # "radius" is used for the algorithms flux.sinc and flux.naive
        calibAlg = measurementConfig.slots.calibFlux
        radius = measurementConfig.algorithms[calibAlg].radius
        apertures = np.array(measurementConfig.algorithms["flux.aperture"].radii)
        if len(np.where(apertures == radius)[0]) == 0:
            raise RuntimeError(
                "Calibration aperture (algorithm %s, radius %f) is not measured by flux.aperture (radii %s)" %
                (calibAlg, radius, apertures))
        calibIndex = np.where(apertures == radius)[0][0]
        corrIndex = np.where(np.isfinite(self.apertureFlux))[0][-1] # Biggest aperture with good correctn

        # The 'ratio' should be less than unity, because we're getting additional flux between
        # the small and large apertures. We will divide fluxes by this ratio so they get brighter.
        # We don't worry about the error except to log it, as it is redundant with any photometric calibration.
        ratio, ratioErr = self.getRatio(calibIndex, corrIndex)
        if log:
            log.info("Applying curve of growth (radius %.1f --> %.1f): %f (+/- %f)" %
                     (apertures[calibIndex], apertures[corrIndex], ratio, ratioErr))

        if algorithms is None:
            algorithms = getApCorrRegistry()
        if catalog is not None:
            for alg in algorithms:
                if alg in catalog:
                    # Fluxes should be divided by the ratio, to get brighter.
                    catalog[alg][:] /= ratio
            metadata = catalog.getMetadata()
            metadata.set("CURVE_OF_GROWTH.CORRECTION.VALUE", ratio)
            metadata.set("CURVE_OF_GROWTH.CORRECTION.ERROR", ratioErr)
            metadata.set("CURVE_OF_GROWTH.CORRECTION.RADIUS.FROM", radius)
            metadata.set("CURVE_OF_GROWTH.CORRECTION.RADIUS.TO", apertures[corrIndex])
            metadata.set("CURVE_OF_GROWTH.CORRECTION.ALGORITHMS", list(algorithms))

        if matches is not None:
            def correctMatches(srcList):
                schema = None
                for src in srcList:
                    if schema is None:
                        schema = src.schema
                    else:
                        assert schema == src.schema, "Schema mismatch"
                for alg in algorithms:
                    if alg in schema:
                        key = schema[alg].asKey()
                        for src in srcList:
                            src[key] /= ratio
            if matchesType is None or matchesType in (self.MATCH_FIRST, self.MATCH_BOTH):
                correctMatches([match.first for match in matches])
            if matchesType is None or matchesType in (self.MATCH_SECOND, self.MATCH_BOTH):
                correctMatches([match.second for match in matches])

        if calib is not None:
            # With this correction our flux measurements are brighter, so to have the same calibrated flux for
            # a reference star as we did previously the fluxMag0 also needs to be larger, so we divide.
            calib /= ratio

        if apCorr is not None:
            for alg in algorithms:
                if alg in apCorr:
                    # Aperture corrections should be divided by the ratio, to get brighter.
                    apCorr[alg] = apCorr[alg]/ratio

class CurveOfGrowth(object):
    def __init__(self, curveOfGrowthCandidateKey=None, curveOfGrowthUsedKey=None,
                 fracInterpolatedMax=0.1, minAnnularFlux=0, nAperture=None, skyNoiseFloor=0.0):
        """Represent a curve of growth generated from a number of stars

        \param curveOfGrowthCandidateKey Key to use to flag candidate objects to use for CurveOfGrowth
        \param curveOfGrowthUsedKey  Key to use to flag objects actually used for CurveOfGrowth
        \param fracInterpolatedMax Maximum fraction of interpolated pixels allowable in an annulus
        \param minAnnularFlux Minimum acceptable flux-per-pixel
        \param nAperture Maximum number of aperture fluxes to use (all, if None)
        \param skyNoiseFloor The per-pixel std. dev. of our knowledge of the sky, to be added in quadrature

        We need to generate curves of growth of stars, and unfortunately stars bright enough to
        measure halo properties are saturated, so a simple addition doesn't suffice. Usage:

           cog = CurveOfGrowth()
           for source in sources:  # for each star that you want included
              cog.addSource(source)

           cog.estimate()
        to estimate a curve of growth from the stars that you added

        Once the curve of growth is estimated it may be returned using
           cog.get("flux.aperture")
           cog.get("flux.aperture.err")
           cog.get("flux.aperture.flags")
           cog.get("flux.aperture.radii")
        A KeyError will be raised if you ask for anything else
        """
        self.curveOfGrowthCandidateKey = curveOfGrowthCandidateKey
        self.curveOfGrowthUsedKey = curveOfGrowthUsedKey
        self.fracInterpolatedMax = fracInterpolatedMax
        self.minAnnularFlux = minAnnularFlux

        self.n = 0
        self.nAperture = nAperture
        self.radii = None
        self.area = None
        self.annularFlux = None
        self.annularFluxErr = None
        self.apertureFlux = None
        self.apertureFluxErr = None
        self.psfFlux = None
        
        self.skyNoiseFloor = skyNoiseFloor
        self.profs = []                 # good profiles, added by addSource
        self.badProfs = []
        self._result = None     # cache of CurveOfGrowthResult object

    @property
    def result(self):
        assert self.apertureFlux is not None and self.apertureFluxErr is not None, \
               "result being used before it has been calculated"
        if self._result is None:
            self._result = CurveOfGrowthResult(self.apertureFlux, self.apertureFluxErr)
        return self._result

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def addSource(self, source):
        """!Add the aperture fluxes in Source to the set to be used in estimating the curve of growth
        \param source -- The source

        We shall ignore radial bins where the Source's aperture flux is missing or too faint to be
        interesting,

        Throw IndexError if the source is useless (<= 1 points)

        """
        # suitably scaled. First make a copy of the input profile/standard deviation

        if self.radii is None:
            md = source.getTable().getMetadata()
            name = "flux_aperture_radii"
            if not md.exists(name):
                raise RuntimeError("You must provide %s in the catalogue metadata" % name)

            self.radii = list(md.get(name)) # outer radius of apertures
            if self.nAperture:
                self.radii = self.radii[0:self.nAperture]

            self.area = np.pi*(np.array(self.radii + [0])**2 -
                               np.array([0] + self.radii)**2)[0:-1] # areas of annuli
            self.radii = np.array(self.radii)

        try:
            self.profs.append(SingleProfile(source, self.fracInterpolatedMax, self.minAnnularFlux,
                                            self.nAperture, self.skyNoiseFloor))
        except IndexError:             # no valid points
            pass

        if self.curveOfGrowthCandidateKey:
            source.set(self.curveOfGrowthCandidateKey, True)

    def _estimate(self, doErrors=True):
        """
        Given a CurveOfGrowth filled with sources' aperture fluxes, derive the curve of growth
        by a straight solution of the MLE problem (i.e. no clipping)
        """
        nProf = len(self.profs)         # number of objects with measured aperture fluxes
        assert nProf > 0, "There must be at least one source to process"

        # determine the radial extent of the data, and the first good value
        i0 = min([p.i0 for p in self.profs])
        nRadial = max([p.npoint for p in self.profs])

        if i0 != 0:
            raise RuntimeError("Profile doesn't extend to r==0")
        
        #
        # ensure that the apertures overlap sufficiently to be able to estimate a
        # curve of growth. conn[k] == 1 => annulus k is connected to annulus k + 1
        #
        conn = np.zeros(nRadial, dtype=int)        # are all radial points connected by data?
        for prof in self.profs:
            conn[prof.i0: prof.npoint] = 1

        disconnected = np.where(conn == 0)[0]
        if disconnected:
            msg = []
            for i in disconnected:
                msg.append("%d-%d" % (i, i+1))
                raise RuntimeError("%s disconnected: %s" % (
                    ("An annulus is" if len(disconnected) == 1 else "Some annuli are"), ", ".join(msg)))

        # w are the inverse-variances
        w = np.empty((nProf, nRadial))
        for j, prof in enumerate(self.profs):
            goodSlice = slice(prof.i0, prof.npoint)

            w[j][0: prof.i0] = 0
            w[j][goodSlice] = (prof.annularFlux[goodSlice]/prof.annularFluxErr[goodSlice])**2

            w[j][prof.npoint:] = 0
        #
        # Setup matrix and vector to solve for nRadial - 1 radial values and nProf scalings alpha.
        #
        # The equation to solve may be written as y = M theta + epsilon
        #
        # ("nRadial - 1" because you can't solve for both the central value of the curve of growth and alpha)
        #
        nN = (nRadial - 1) + nProf

        MTVM = np.empty((nN, nN))          # normal equations matrix (M^T V^-1 M)
        MVy = np.empty(nN)                 # vector M^T V^-1 y
        #
        # In what follows, i is an index over radial bins, and j an index of input sources
        # Remember that we are _not_ estimating self.annularFlux[0]
        #
        for i in range(1, nRadial):     # top left block
            MTVM[i - 1, 0:nRadial - 1] = 0.0
            MTVM[i - 1, i - 1] = np.sum(w[:, i])

        for j in range(nProf):          # bottom right block
            MTVM[nRadial - 1 + j, nRadial - 1: nRadial - 1 + nProf] = 0
            MTVM[nRadial - 1 + j, nRadial - 1 + j] = np.sum(w[j, :])

        for i in range(1, nRadial):      # upper right and lower left blocks
            for j in range(nProf):
                MTVM[nRadial - 1 + j, i - 1] = w[j, i]
                MTVM[i - 1, nRadial - 1 + j] = w[j, i]

        for i in range(1, nRadial):     # top of vector (log(P_i))
            s = 0
            for j, prof in enumerate(self.profs):
                if w[j, i] > 0:
                    s += np.log(abs(prof.annularFlux[i]))*w[j, i]
            MVy[i - 1] = s

        for j, prof in enumerate(self.profs): # bottom of vector (log(alpha^j))
            s = 0
            for i in range(nRadial):
                if w[j, i] > 0:
                    s += np.log(abs(prof.annularFlux[i]))*w[j, i]

            MVy[nRadial - 1 + j] = s
        #
        # OK, invert the matrix and then solve for our estimates. We need the
        # inverse for error estimates, so we don't solve the matrix equation directly
        #
        try:
            if doErrors:
                iMTVM = np.linalg.inv(MTVM)
                theta = np.exp(iMTVM.dot(MVy))
                thetaErr = np.sqrt(np.diag(iMTVM))*theta
            else:
                theta = np.exp(np.linalg.solve(MTVM, MVy))
                thetaErr = np.ones_like(theta)
        except Exception as e:
            raise RuntimeError("Curve of growth matrix is singular: %s" % e)

        self.annularFlux = np.empty(nRadial)
        self.annularFluxErr = np.empty_like(self.annularFlux)

        self.annularFlux[0] = 1
        self.annularFluxErr[0] = 0

        self.annularFlux[1:] = theta[0:nRadial-1]
        self.annularFluxErr[1:] = thetaErr[0:nRadial-1]
        alpha = theta[nRadial - 1:]
        alphaErr = thetaErr[nRadial - 1:]

        for prof, a, aErr in zip(self.profs, alpha, alphaErr):
            prof.alpha, prof.alphaErr = a, aErr
        #
        # Renormalise the cumulative profile to reach 1.0
        #
        fluxMax = np.cumsum(self.annularFlux)[-1]
        self.annularFlux /= fluxMax
        self.annularFluxErr /= fluxMax

        for prof in self.profs:
            prof.alpha *= fluxMax
            prof.alphaErr *= fluxMax
        #
        # set aperture.flux from the annularFlux
        #
        self.apertureFlux = np.cumsum(self.annularFlux)
        self.apertureFluxErr = np.sqrt(np.cumsum(self.annularFluxErr**2))
        #
        # estimate the PSF counts corresponding to a profile with alpha == 1
        #
        self.psfFlux = 0
        s = 0
        for prof in self.profs:
            if prof.i0 == 0:            # not saturated
                s += prof.alpha
                self.psfFlux += prof.psfFlux

        if s == 0:
            self.psfFlux = None
        else:
            self.psfFlux /= s

    def estimate(self, maxRChi2=[100], nSaturated=None, nNonSaturated=None, finalEstimationAlgorithm="mean"):
        """!Given a CurveOfGrowth filled with aperture photometry from sources, derive the curve of growth

        \param maxRChi2 clip the input objects with a reduced chi^2 larger than this (may be a list)
        \param nSaturated How many of the brightest saturated candidates to use
        \param nNonSaturated How many of the brightest non-saturated candidates to use
        \param finalEstimationAlgorithm How to estimate the final apertureFluxes (None, "mean", "median")

        The final returned apertureFlux may be the one estimated from MLE fit, or the bin-by-bin
        (weighted) mean or median; the choice is set by finalEstimationAlgorithm.  Empirically, better results
        seem to be returned by using mean or median.
        """
        #
        # Choose the nSaturated brightest saturated and the brightest nNonSaturated non-saturated candidates
        #
        if nSaturated or nNonSaturated:
            try:
                saturatedKey = self.profs[0].source.getSchema().find("flags.pixel.saturated.center").getKey()
            except KeyError:
                saturatedKey = None

            saturated, nonSaturated = [], []
            for prof in self.profs:
                s = prof.source
                if saturatedKey and s.get(saturatedKey):
                    saturated.append([prof, s.getApFlux()])
                else:
                    nonSaturated.append([prof, s.getPsfFlux()])

            def byFlux(a, b):
                return cmp(a[1], b[1])

            saturated.sort(byFlux)
            if nSaturated:
                saturated = saturated[:nSaturated]

            nonSaturated.sort(byFlux)
            if nNonSaturated:
                nonSaturated = nonSaturated[:nNonSaturated]

            self.profs = [prof for prof, flux in (saturated + nonSaturated)]

        try:
            maxRChi2[0]
        except IndexError:
            maxRChi2 = [maxRChi2]

        #
        # Make an initial estimate without clipping
        #
        self._estimate(doErrors=(len(maxRChi2) == 0))

        #
        # Now estimate with clipping
        #
        for i, maxRChi2Value in enumerate(maxRChi2):
            #
            # Estimate a robust average of the input aperture fluxes in each annulus
            #
            robustAnnularFlux = self._estimateMedianAnnularFlux()
            #
            goodProfs = []
            for prof in self.profs:
                goodSlice = slice(prof.i0, prof.npoint)

                annErr = np.hypot(prof.annularFluxErr[goodSlice], prof.alphaErr*robustAnnularFlux[goodSlice])

                chi = (prof.annularFlux[goodSlice] - prof.alpha*robustAnnularFlux[goodSlice])/annErr
                chi2 = np.sum(chi**2)

                self.rchi2 = chi2/(prof.npoint - prof.i0 - 1)

                if self.rchi2 > maxRChi2Value or not np.all(np.isfinite(chi)):
                    self.badProfs.append(prof)
                else:
                    goodProfs.append(prof)

            self.profs = goodProfs
            #
            # Estimate the curve of growth having discarded those suspect objects
            #
            self._estimate(doErrors=(i + 1 == len(maxRChi2)))
        #
        # Record which sources were actually used
        #
        if self.curveOfGrowthUsedKey:
            for prof in self.profs:
                prof.source.set(self.curveOfGrowthUsedKey, True)
            for prof in self.badProfs:
                prof.source.set(self.curveOfGrowthUsedKey, False)

        if finalEstimationAlgorithm:
            if finalEstimationAlgorithm == "mean":
                robustAnnularFlux = self._estimateMeanAnnularFlux()
            elif finalEstimationAlgorithm == "median":
                robustAnnularFlux = self._estimateMedianAnnularFlux()
            else:
                raise RuntimeError("Unknown value for finalEstimationAlgorithm: %s" %
                                   finalEstimationAlgorithm)

            scale = np.cumsum(self.annularFlux)[-1]/np.cumsum(robustAnnularFlux)[-1]

            self.annularFlux = scale*robustAnnularFlux
            self.annularFluxErr *= scale
            self.psfFlux *= scale

            for prof in self.profs:
                prof.alpha /= scale

        self.radii = self.radii[0:len(self.annularFlux)] # maybe none of the profiles were complete

    def get(self, what):
        """!Return the estimated aperture flux and its error
        \param what "flux.aperture" or "flux.aperture.err"
        """
        try:
            return {"flux.aperture"       : self.apertureFlux,
                    "flux.aperture.err"   : self.apertureFluxErr,
                    "flux.aperture.flags" : self.apertureFlux is None,
                    "flux.aperture.radii" : self.radii,
            }[what]
        except NameError, e:
            raise KeyError(str(e))      # be consistent with afwTable

    def getRatio(self, inner, outer):
        """!Return the flux ratio between two radii and the error on the ratio.

        \param inner     Index of the inner radius (numerator of the ratio); if -ve relative to end of array
        \param outer     Index of the outer radius (denominator of the ratio); if -ve relative to end of array
        """
        return self.result.getRatio(inner, outer)

    def _estimateMeanAnnularFlux(self):
        """
        Given a CurveOfGrowth filled with aperture photometry from sources, return the weighted mean
        of the input aperture fluxes scaled by the estimated alpha values
        """
        meanAnnularFlux = np.zeros_like(self.annularFlux)
        for i in range(len(meanAnnularFlux)):
            values, valuesErr = [], []
            for prof in self.profs:
                if i >= prof.i0 and i < prof.npoint: # a valid point
                    values.append(prof.annularFlux[i]/prof.alpha)
                    valuesErr.append(prof.annularFluxErr[i]/prof.alpha)

            meanAnnularFlux[i] = np.average(values, weights=1/np.array(valuesErr)**2)

        return meanAnnularFlux

    def _estimateMedianAnnularFlux(self):
        """
        Given a CurveOfGrowth filled with aperture photometry from sources, return the median
        of the input aperture fluxes scaled by the estimated alpha values
        """
        medianAnnularFlux = np.zeros_like(self.annularFlux)
        for i in range(len(medianAnnularFlux)):
            values = []
            for prof in self.profs:
                if i >= prof.i0 and i < prof.npoint: # a valid point
                    values.append(prof.annularFlux[i]/prof.alpha)

            medianAnnularFlux[i] = np.median(values)

        return medianAnnularFlux

    def estimatePsfFluxFromApertureFlux(self,
                                        source,     # the object of interest
                                        i0=0,       # and the first good radial point
                                        rMax=None   # use only annuli up to this radius
    ):
        """
        Calculate the psfFlux given a set of aperture fluxes for a source.

        This is done by matching the profiles (so the source can perfectly well be saturated).
        
        We could use a technique similar to that used to derive the curve of growth,
        but let us be lazy and assume that the errors in the curve of growth are small
        """

        nRadial = source.get("flux.aperture.nProfile")
        if nRadial == 0:
            raise IndexError("Source %d has no aperture measurements" % source.getId())
        
        flux = source.get("flux.aperture")
        fluxErr = source.get("flux.aperture.err")

        if len(self.annularFlux) < nRadial:
            nRadial = len(self.annularFlux)

        if rMax is not None:            # maximum radius to use
            for i in range(i0, nRadial):
                if self.radii[i + 1] > rMax:
                    break

                if i > i0:                      # don't entirely trim overlap
                    nRadial = i

        if i0 >= nRadial:
            raise RuntimeError("Available annular fluxes and CurveOfGrowth have no points in common")

        sumComp = 0
        sumObj = 0

        prevFlux = flux[i0 - 1] if i0 > 1 else 0.0
        for i in range(i0, nRadial):
            iannularFluxErr2 = 1/(fluxErr[i]**2 + self.annularFluxErr[i]**2) # inverse variance
            sumObj  += (flux[i] - prevFlux)*self.annularFlux[i]   *iannularFluxErr2
            sumComp +=                      self.annularFlux[i]**2*iannularFluxErr2

            prevFlux = flux[i]
            
        assert sumComp != 0, "The curve of growth must have net flux"

        return self.psfFlux*sumObj/sumComp

    def plot(self, normalize=True, showRadialProfile=False, alpha=0.1, yscaleType=None,
             minLabelledValue=None, showInputs=True, clf=True, title=None):
        """!Plot the fit curve of growth and the input aperture fluxes using matplotlib
        
        \param normalize If True, normalize all the constituent sources to constant flux
        \param showRadialProfile If True, plot the surface brightness rather than the curve of growth
        \param alpha transparency for matplotlib
        \param yscaleType How to plot y axis (None, 'log' or 'linear')
        \param minLabelledValue Curves with final values larger than minLabelledValue are labelled
        \param showInputs Show the objects used to construct the curve of growth
        \param clf Clear the figure before plotting
        \param title A title to use (None: number of input objects)
        """
        mpAlpha = alpha                 # we use the name alpha for other things

        if not plt:
            print >> sys.stderr, "Failed to import matplotlib.pyplot"
            return

        global fig
        if not fig:
            fig = plt.figure()
        elif clf:
            fig.clf()

        axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

        if showRadialProfile:
            if not yscaleType:
                yscaleType = 'log'
            def profileFunction(annularFlux, annularFluxErr=None, prof=None, cogAnnularFlux=None):
                profile = annularFlux/self.area

                if annularFluxErr is None:
                    return profile
                else:
                    return profile, annularFluxErr/self.area
        else:
            if not yscaleType:
                yscaleType = 'linear'

            def profileFunction(annularFlux, annularFluxErr=None, prof=None, cogAnnularFlux=None):
                if prof and prof.i0 > 0 and cogAnnularFlux is not None:
                    annularFlux = annularFlux[:] 
                    annularFlux[0:prof.i0] = prof.alpha*cogAnnularFlux[0:prof.i0]

                profile = np.cumsum(annularFlux)

                if annularFluxErr is None:
                    return profile
                else:
                    return profile, np.cumsum(annularFluxErr)

        robustFlux = profileFunction(self._estimateMedianAnnularFlux())

        if normalize:
            flux, fluxErr = profileFunction(self.annularFlux, self.annularFluxErr)
            scale = 1.0/flux[0] if showRadialProfile else 1.0

            robustFlux *= scale
            flux *= scale
            fluxErr *= scale
            
            axes.plot(self.radii, flux, '-', color='black')

            nr = len(self.radii)
            r = np.empty(2*nr)
            r[:nr] = self.radii
            r[nr:] = list(reversed(self.radii))

            fluxPM = np.empty_like(r)
            fluxPM[:nr] = flux + fluxErr
            fluxPM[nr:] = list(reversed(flux - fluxErr))

            if yscaleType == 'log':     # we shouldn't have to do this...
                fluxPM = np.where(fluxPM > 0, np.log(fluxPM), np.log(1e-6))

            p = plt.Polygon(zip(r, fluxPM), alpha=0.4)
            axes.add_artist(p)

            if not showRadialProfile:
                axes.axhline(1.0, ls=':')

            axes.errorbar(self.radii, robustFlux, yerr=fluxErr, ls=':', color='black')

        nlabel = 0
        for i, prof in enumerate(sorted(self.profs, lambda a, b: cmp(b.alpha, a.alpha)), -len(self.profs)/2):
            if not showInputs:
                continue

            good = slice(prof.i0, prof.npoint)

            flux, fluxErr = profileFunction(prof.annularFlux, prof.annularFluxErr, prof, self.annularFlux)
            if normalize:
                flux *= scale/prof.alpha
                fluxErr *= scale/prof.alpha

            label = None
            if minLabelledValue is not None and flux[good][-1] > minLabelledValue:
                label = prof.source.getId() & 0xffff
            lines = axes.errorbar(self.radii[good] + 0.1/len(self.profs)*i, flux[good], yerr=fluxErr[good],
                                  fmt='o', label=label, alpha=mpAlpha)
            nlabel += (label != None)

            color = lines[0].get_color()
            axes.plot(self.radii[good], flux[good], ls='-', color=color, alpha=mpAlpha)

            if not normalize:
                p, pErr = profileFunction(self.annularFlux, self.annularFluxErr)

                axes.errorbar(self.radii[good] + 0.0,
                             prof.alpha*p[good], yerr=prof.alpha*pErr[good],
                             ls='--', color=color)

                axes.errorbar(self.radii[good],
                             prof.alpha*robustFlux[good], yerr=fluxErr[good],
                             ls=':', color=color)

        if 0 < nlabel < 10:
            legend = axes.legend(loc='best')
            if legend:
                legend.draggable(True)

        axes.set_xlabel("Radius (pixels)")
        if showRadialProfile:
            ylabel = "Radial profile"
        else:
            ylabel = "Cumulative flux"

        if normalize:
            ylabel += " (normalised)"

        axes.set_ylabel(ylabel)
        axes.set_yscale(yscaleType)
        if title is None:
            title = "%d objects" % len(self.profs)
        axes.set_title(title)

        return fig
    
    def apply(self, *args, **kwargs):
        """!Apply the results of the curve of growth analysis to a source catalog

        We deduce the correction that needs to be applied from inspecting the
        'measurementConfig' (a SourceMeasurementConfig) and finding the aperture
        being used for calibration.  The correction is then applied to the flux
        and error of each of the specified 'algorithms'. If 'algorithms' is None,
        the set of algorithms in the aperture correction registry is used.

        \param measurementConfig    Configuration for source measurement
        \param catalog    Source catalog to have flux measurements corrected
        \param algorithms    Iterable of algorithms to correct, or None
        \param calib    Calib object to be corrected
        \param apCorr    Aperture corrections to be corrected
        \param log    Log object for logging the factor that's applied
        """
        return self.result.apply(*args, **kwargs)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SingleProfile(object):
    def __init__(self, source, fracInterpolatedMax, minAnnularFlux, nAperture=None, skyNoiseFloor=0):
        """!The aperture and annular fluxes of a single Source, which define points on a radial profile

        \param source  A measured Source
        \param fracInterpolatedMax Maximum fraction of interpolated pixels allowable in an annulus
        \param minAnnularFlux Minimum acceptable flux-per-pixel
        \param nAperture Maximum number of aperture fluxes to use (all, if None)
        \param skyNoiseFloor  The per-pixel std. dev. of our knowledge of the sky, to be added in quadrature
        """
        radii = source.getTable().getMetadata().get("flux_aperture_radii")
        #
        # We assume that the radii are monotonic increasing, so let's check
        #
        assert len(radii), "You must provide at least one radial point"
        r0 = radii[0]
        for r in radii[1:]:
            assert r > r0, "aperture radii must be monotonic increasing"
            r0 = r

        flux = source.get("flux.aperture")
        fluxErr = source.get("flux.aperture.err")
        nInterpPixel = source.get("flux.aperture.nInterpolatedPixel")

        if nAperture:
            radii = radii[0:nAperture]
            flux = flux[0:nAperture]
            fluxErr = fluxErr[0:nAperture]
            nInterpPixel = nInterpPixel[0:nAperture]

        self.i0 = None                  # first good value
        self.annularFlux = np.empty_like(radii) + np.nan    # star's mean annular fluxes
        self.annularFluxErr = np.empty_like(radii) # errors in annularFlux
        self.annularArea = np.empty_like(radii)    # areas of annuli

        areaOld = 0                   # area of previous aperture
        fluxOld = 0                   # flux inside inside areaOld
        fluxErrOld = 0                # error in flux inside inside areaOld
        nInterpOld = 0                # number of interpolated pixels inside areaOld

        nRadial = source.get("flux.aperture.nProfile")
        if nAperture and nRadial > nAperture:
            nRadial = nAperture

        if nRadial == 0:
            raise IndexError("Source %d has no aperture measurements" % source.getId())
        
        self.prof = flux[:]

        self.npoint = 0         # (number of good values) + i0
        for i in range(nRadial):
            if not np.isfinite(flux[i]) or not np.isfinite(fluxErr[i]) or fluxErr[i] < fluxErrOld:
                break
            self.annularFlux[i] = flux[i] - fluxOld
            self.annularFluxErr[i] = np.sqrt(fluxErr[i]**2 - fluxErrOld**2)

            nInterpAperture = nInterpPixel[i] - nInterpOld
            area = np.pi*radii[i]**2
            annularArea = area - areaOld
            self.annularArea[i] = annularArea

            fracInteropolated = nInterpAperture/annularArea
            if fracInteropolated < fracInterpolatedMax: # acceptable annulus
                if self.i0 is None:
                    self.i0 = i
            
            if self.annularFlux[i] < minAnnularFlux*annularArea:
                break

            areaOld, fluxOld, fluxErrOld, nInterpOld = area, flux[i], fluxErr[i], nInterpPixel[i]
            self.npoint += 1

        self.source = source             # debugging
        self.psfFlux = source.getPsfFlux()
        self.alpha, self.alphaErr = None, None

        if self.i0 is None or self.npoint <= self.i0 + 1:   # no new information (1 point's useless)
            raise IndexError("Source %d has <= 1 valid aperture measurement" % source.getId())
        #
        # Impose a minimum variance due to unmodelled sky noise
        #
        if skyNoiseFloor:
            goodSlice = slice(self.i0, self.npoint)
            self.annularFluxErr[goodSlice] = np.sqrt(self.annularFluxErr[goodSlice]**2 + 
                                                     self.annularArea[goodSlice]*skyNoiseFloor**2)

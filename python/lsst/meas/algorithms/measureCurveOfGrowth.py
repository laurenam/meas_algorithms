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

try:
    import matplotlib.pyplot as plt
    fig = None
except ImportError as e:
    print e
    plt = None

__all__ = ("CurveOfGrowthMeasurementConfig", "CurveOfGrowthMeasurementTask", "CurveOfGrowth")

class CurveOfGrowthMeasurementConfig(pex_config.Config):
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
        default = 0.5,
        check = lambda x: x >= 0.0,
    )
    psfFluxMin = pex_config.Field(
        doc = "Minimum value of psf.flux for an object to be included in curve of growth",
        dtype = float,
        default = 5e4,
        check = lambda x: x >= 0.0,
    )
    maxRChi2 = pex_config.ListField(
        doc = """List of values of reduced chi^2 that should be applied in order to clip sources

        N.b. these values are so large because of contamination in the annuli by to e.g. the faint wings of
        neighbouring objects.  The real solution here is to write cleverer aperture flux code, \'a la SDSS
        """,
        dtype = float,
        default = [100, 75, 50, 20],
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

    def run(self, catalog):
        """!Estimate the curve of growth from a set of Sources by invoking run(catalog)
        \param catalog  A catalog of Sources

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
        sch = catalog.getSchema()

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
                            self.config.fracInterpolatedMax, self.config.minAnnularFlux)

        for s in catalog:
            if (deblend_nchild is not None and s.get(deblend_nchild) > 0) or \
               s.getPsfFlux() < self.config.psfFluxMin or \
               s.get(classification_extendedness) > self.config.classificationMax:
                continue

            for badFlag in badFlags:
                if s.get(badFlag):
                    continue

            cog.addSource(s)
        #
        # and actually estimate it
        #
        cog.estimate(maxRChi2=self.config.maxRChi2,
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

class CurveOfGrowth(object):
    def __init__(self, curveOfGrowthCandidateKey=None, curveOfGrowthUsedKey=None,
                 fracInterpolatedMax=0.1, minAnnularFlux=0):
        """Represent a curve of growth generated from a number of stars

        \param curveOfGrowthCandidateKey Key to use to flag candidate objects to use for CurveOfGrowth
        \param curveOfGrowthUsedKey  Key to use to flag objects actually used for CurveOfGrowth
        \param fracInterpolatedMax Maximum fraction of interpolated pixels allowable in an annulus
        \param minAnnularFlux Minimum acceptable flux-per-pixel

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
        self.radii = None
        self.area = None
        self.annularFlux = None
        self.annularFluxErr = None
        self.apertureFlux = None
        self.apertureFluxErr = None
        self.psfFlux = None
        
        self.profs = []                 # good profiles, added by addSource
        self.badProfs = []

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
            if md.exists(name):
                self.radii = list(md.get(name)) # outer radius
            else:
                raise RuntimeError("You must provide %s in the catalogue metadata" % name)

            self.area = np.pi*(np.array(self.radii + [0])**2 -
                               np.array([0] + self.radii)**2)[0:-1] # areas of annuli
            self.radii = np.array(self.radii)

        try:
            self.profs.append(SingleProfile(source, self.fracInterpolatedMax, self.minAnnularFlux))
        except IndexError:             # no valid points
            pass

        if self.curveOfGrowthCandidateKey:
            source.set(self.curveOfGrowthCandidateKey, True)

    def _estimate(self):
        """
        Given a CurveOfGrowth filled with sources' aperture fluxes, derive the curve of growth
        by a straight solution of the MLE problem (i.e. no clipping)
        """
        nProf = len(self.profs)         # number of objects with measured aperture fluxes
        assert nProf > 0, "There is at least one source to process"

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
            iMTVM = np.linalg.inv(MTVM)
        except Exception as e:
            raise RuntimeError("Curve of growth matrix is singular: %s" % e)

        self.annularFlux = np.empty(nRadial)
        self.annularFluxErr = np.empty_like(self.annularFlux)

        self.annularFlux[0] = 1
        self.annularFluxErr[0] = 0

        theta = np.exp(iMTVM.dot(MVy))
        thetaErr = np.sqrt(np.diag(iMTVM))*theta

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

    def estimate(self, maxRChi2=[100], finalEstimationAlgorithm="mean"):
        """!Given a CurveOfGrowth filled with aperture photometry from sources, derive the curve of growth

        \param maxRChi2 clip the input objects with a reduced chi^2 larger than this (may be a list)
        \param finalEstimationAlgorithm How to estimate the final apertureFluxes (None, "mean", "median")

        The final returned apertureFlux may be the one estimated from MLE fit, or the bin-by-bin
        (weighted) mean or median; the choice is set by finalEstimationAlgorithm.  Empirically, better results
        seem to be returned by using mean or median.
        """
        #
        # Make an initial estimate without clipping
        #
        self._estimate()

        try:
            maxRChi2[0]
        except IndexError:
            maxRChi2 = [maxRChi2]

        for maxRChi2 in maxRChi2:
            #
            # Estimate a robust average of the input aperture fluxes in each annulus
            #
            robustAnnularFlux = self._estimateMedianAnnularFlux()
            #
            goodProfs = []
            for prof in self.profs:
                goodSlice = slice(prof.i0, prof.npoint)
                chi = ((prof.annularFlux[goodSlice] - prof.alpha*robustAnnularFlux[goodSlice])/
                       prof.annularFluxErr[goodSlice])
                chi2 = np.sum(chi**2)

                self.rchi2 = chi2/(prof.npoint - prof.i0 - 1)

                if self.rchi2 > maxRChi2 or not np.all(np.isfinite(chi)):
                    self.badProfs.append(prof)
                else:
                    goodProfs.append(prof)

            self.profs = goodProfs
            #
            # Estimate the curve of growth having discarded those suspect objects
            #
            self._estimate()
        #
        # Record which sources were actually used
        #
        if self.curveOfGrowthUsedKey:
            for prof in self.profs:
                prof.source.set(self.curveOfGrowthUsedKey, True)

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

        \param inner     Index of the inner radius (numerator of the ratio).
        \param outer     Index of the outer radius (denominator of the ratio).
        """
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

    def plot(self, normalize=True, showRadialProfile=False, alpha=0.1):
        """!Plot the fit curve of growth and the input aperture fluxes using matplotlib
        
        \param normalize If True, normalize all the constituent sources to constant flux
        \param showRadialProfile If True, plot the surface brightness rather than the curve of growth
        \param alpha transparency for matplotlib
        """
        mpAlpha = alpha                 # we use the name alpha for other things

        if not plt:
            print >> sys.stderr, "Failed to import matplotlib.pyplot"
            return

        global fig
        if not fig:
            fig = plt.figure()
        else:
            fig.clf()

        axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

        if showRadialProfile:
            def profileFunction(annularFlux, annularFluxErr=None):
                prof = annularFlux/self.area

                if annularFluxErr is None:
                    return prof
                else:
                    return prof, annularFluxErr/self.area
        else:
            def profileFunction(annularFlux, annularFluxErr=None):
                prof = np.cumsum(annularFlux)

                if annularFluxErr is None:
                    return prof
                else:
                    return prof, np.cumsum(annularFluxErr)

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

            p = plt.Polygon(zip(r, fluxPM), alpha=0.4)
            axes.add_artist(p)

            if not showRadialProfile:
                axes.axhline(1.0, ls=':')

            axes.errorbar(self.radii, robustFlux, yerr=fluxErr, ls=':', color='black')

        nlabel = 0
        for i, prof in enumerate(sorted(self.profs, lambda a, b: cmp(b.alpha, a.alpha)), -len(self.profs)/2):
            good = slice(prof.i0, prof.npoint)

            flux, fluxErr = profileFunction(prof.annularFlux, prof.annularFluxErr)
            if normalize:
                flux *= scale/prof.alpha
                fluxErr *= scale/prof.alpha

            label = prof.source.getId() if (flux[good][-1] > (10 if normalize else 4e5)) else None
            lines = axes.errorbar(self.radii[good]+0.01*i, flux[good], yerr=fluxErr[good],
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

            axes.set_xlim(None, 75)
            
        axes.set_xlabel("Radius (pixels)")
        if showRadialProfile:
            ylabel = "Radial profile"
        else:
            ylabel = "Cumulative flux"

        if normalize:
            ylabel += " (normalised)"

        axes.set_ylabel(ylabel)

        return fig
    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SingleProfile(object):
    def __init__(self, source, fracInterpolatedMax, minAnnularFlux):
        """!The aperture and annular fluxes of a single Source, which define points on a radial profile

        \param source  A measured Source
        \param fracInterpolatedMax Maximum fraction of interpolated pixels allowable in an annulus
        \param minAnnularFlux Minimum acceptable flux-per-pixel
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

        self.i0 = None                  # first good value
        self.annularFlux = np.empty_like(radii) + np.nan    # star's mean annular fluxes
        self.annularFluxErr = np.empty_like(radii) # errors in annularFlux
        self.annularArea = np.empty_like(radii)    # areas of annuli

        areaOld = 0                   # area of previous aperture
        fluxOld = 0                   # flux inside inside areaOld
        fluxErrOld = 0                # error in flux inside inside areaOld
        nInterpOld = 0                # number of interpolated pixels inside areaOld

        nRadial = source.get("flux.aperture.nProfile")
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

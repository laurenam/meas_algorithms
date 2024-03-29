# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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
import math
import numpy as np
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
import lsst.afw.display.ds9 as ds9

from . import algorithmsLib
from .algorithmRegistry import *
from .replaceWithNoise import *

__all__ = "SourceSlotConfig", "SourceMeasurementConfig", "SourceMeasurementTask"

class SourceSlotConfig(pexConfig.Config):

    centroid = pexConfig.Field(dtype=str, default="centroid.sdss", optional=True,
                             doc="the name of the centroiding algorithm used to set source x,y")
    shape = pexConfig.Field(dtype=str, default="shape.sdss", optional=True,
                          doc="the name of the algorithm used to set source moments parameters")
    apFlux = pexConfig.Field(dtype=str, default="flux.sinc", optional=True,
                           doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = pexConfig.Field(dtype=str, default="flux.gaussian", optional=True,
                           doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = pexConfig.Field(dtype=str, default="flux.psf", optional=True,
                            doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = pexConfig.Field(dtype=str, default="flux.gaussian", optional=True,
                             doc="the name of the algorithm used to set the source inst flux slot")
    calibFlux = pexConfig.Field(dtype=str, default="flux.naive", optional=True,
                                doc="the name of the aperture flux algorithm used for calibration")

    def setupTable(self, table, prefix=None):
        """Convenience method to setup a table's slots according to the config definition.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        if prefix is None: prefix = ""
        if self.centroid is not None: table.defineCentroid(prefix + self.centroid)
        if self.shape is not None: table.defineShape(prefix + self.shape)
        if self.apFlux is not None: table.defineApFlux(prefix + self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(prefix + self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(prefix + self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(prefix + self.instFlux)
        if self.calibFlux is not None: table.defineCalibFlux(prefix + self.calibFlux)

class SourceMeasurementConfig(pexConfig.Config):
    """
    Configuration for SourceMeasurementTask.
    A configured instance of MeasureSources can be created using the
    makeMeasureSources method.
    """

    slots = pexConfig.ConfigField(
        dtype = SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source.\n"
        )

    algorithms = AlgorithmRegistry.all.makeField(
        multi=True,
        default=["flags.pixel",
                 "centroid.naive",
                 "shape.sdss",
                 "flux.gaussian", "flux.naive", "flux.psf", "flux.sinc",
                 "correctfluxes",
                 "classification.extendedness",
                 "skycoord",
                 ],
        doc="Algorithms that will be run by default."
        )

    centroider = AlgorithmRegistry.filter(CentroidConfig).makeField(
        multi=False, default="centroid.sdss", optional=True,
        doc="Configuration for the initial centroid algorithm used to\n"\
            "feed center points to other algorithms.\n\n"\
            "Note that this is in addition to the centroider listed in\n"\
            "the 'algorithms' field; the same name should not appear in\n"\
            "both.\n\n"\
            "This field DOES NOT set which field name will be used to define\n"\
            "the alias for source.getX(), source.getY(), etc.\n"
        )

    doBlendedness = pexConfig.Field(dtype=bool, default=False,
                                    doc="Whether to compute blendedness metrics")

    blendedness = pexConfig.ConfigField(dtype=algorithmsLib.BlendednessConfig,
                                        doc="Configuration for blendedness metrics")

    doReplaceWithNoise = pexConfig.Field(dtype=bool, default=True, optional=False,
                                         doc='When measuring, replace other detected footprints with noise?')

    replaceWithNoise = pexConfig.ConfigurableField(
        target = ReplaceWithNoiseTask,
        doc = ("Task for replacing other sources by noise when measuring sources; run when " +
               "'doReplaceWithNoise' is set."),
    )

    prefix = pexConfig.Field(dtype=str, optional=True, default=None, doc="prefix for all measurement fields")

    def validate(self):
        pexConfig.Config.validate(self)
        if self.centroider.name in self.algorithms.names:
            raise ValueError("The algorithm in the 'centroider' field must not also appear in the "\
                                 "'algorithms' field.")
        if self.slots.centroid is not None and (self.slots.centroid not in self.algorithms.names
                                                and self.slots.centroid != self.centroider.name):
            raise ValueError("source centroid slot algorithm '%s' is not being run." % self.slots.astrom)
        if self.slots.shape is not None and self.slots.shape not in self.algorithms.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux, self.slots.instFlux):
            if slot is not None:
                for name in self.algorithms.names:
                    if len(name) <= len(slot) and name == slot[:len(name)]:
                        break
                else:
                    raise ValueError("source flux slot algorithm '%s' is not being run." % slot)
                

    def makeMeasureSources(self, schema, metadata=None, isForced=False):
        """ Convenience method to make a MeasureSources instance and
        fill it with the configured algorithms.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        builder = algorithmsLib.MeasureSourcesBuilder(
            self.prefix if self.prefix is not None else "",
            isForced
            )
        if self.centroider.name is not None:
            builder.setCentroider(self.centroider.apply())
        builder.addAlgorithms(self.algorithms.apply())
        return builder.build(schema, metadata)

class SourceMeasurementTask(pipeBase.Task):
    """Measure the properties of sources on a single exposure.

    This task has no return value; it only modifies the SourceCatalog in-place.
    """
    ConfigClass = SourceMeasurementConfig
    _DefaultName = "sourceMeasurement"

    def __init__(self, schema, algMetadata=None, isForced=False, **kwds):
        """Create the task, adding necessary fields to the given schema.

        @param[in,out] schema        Schema object for measurement fields; will be modified in-place.
        @param[in,out] algMetadata   Passed to MeasureSources object to be filled with initialization
                                     metadata by algorithms (e.g. radii for aperture photometry).
        @param[in]     isForced      If true, measurements will be performed in "forced" mode,
                                     and the references and refWcs arguments must be passed to run()
                                     and measure().
        @param         **kwds        Passed to Task.__init__.
        """
        pipeBase.Task.__init__(self, **kwds)
        self.measurer = self.config.makeMeasureSources(schema, algMetadata, isForced=isForced)
        if self.config.doReplaceWithNoise:
            self.makeSubtask('replaceWithNoise')
        if self.config.doBlendedness:
            self.blendedness = algorithmsLib.Blendedness(self.config.blendedness.makeControl(), schema)

    def preMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the start of the
        measure() method.'''

        # pipe_base's Task provides self._display.
        if self._display:
            frame = 0
            ds9.mtv(exposure, title="input", frame=frame)

    def postMeasureHook(self, exposure, sources):
        '''A hook, for debugging purposes, that is called at the end of the
        measure() method.'''

    def preSingleMeasureHook(self, exposure, sources, i):
        '''A hook, for debugging purposes, that is called immediately
        before the measurement algorithms for each source.

        Note that this will also be called with i=-1 just before entering the
        loop over measuring sources.'''

        if self._display:
            if i < 0:
                # First time...
                try:
                    self.deblendAsPsfKey = sources.getSchema().find("deblend.deblended-as-psf").getKey()
                except KeyError:
                    self.deblendAsPsfKey = None
                    
            if self._display > 2 and i >= 0:
                peak = sources[i].getFootprint().getPeaks()[0]
                print sources[i].getId(), peak.getIx(), peak.getIy()

    def postSingleMeasureHook(self, exposure, sources, i):
        '''A hook, for debugging purposes, that is called immediately after
        the measurement algorithms.'''
        self.postSingleMeasurementDisplay(exposure, sources[i])

    def postSingleMeasurementDisplay(self, exposure, source):
        if self._display:
            if self._display > 1:
                ds9.dot(str(source.getId()), source.getX() + 2, source.getY(),
                        size=3, ctype=ds9.RED)
                cov = source.getCentroidErr()
                ds9.dot(("@:%.1f,%.1f,%1f" % (cov[0,0], cov[0,1], cov[0,0])),
                        *source.getCentroid(), size=3, ctype=ds9.RED)
                symb = "%d" % source.getId()
            else:
                symb = "*" if self.deblendAsPsfKey and source.get(self.deblendAsPsfKey) else "+"
                ds9.dot(symb, *source.getCentroid(), size=3,
                        ctype=ds9.RED if source.get("parent") == 0 else ds9.MAGENTA)

                for p in source.getFootprint().getPeaks():
                    ds9.dot("+", *p.getF(), size=0.5, ctype=ds9.YELLOW)
    
    @pipeBase.timeMethod
    def measure(self, exposure, sources, noiseImage=None, noiseMeanVar=None, references=None, refWcs=None,
                beginPriority=0.0, endPriority=float("inf")):
        """Measure sources on an exposure, with no aperture correction.

        @param[in]     exposure       Exposure to process
        @param[in,out] sources        SourceCatalog containing sources detected on this exposure.
        @param[in]     noiseImage     If 'config.doReplaceWithNoise = True', you can pass in
                                      an Image containing noise.  This overrides the "config.noiseSource"
                                      setting.
        @param[in]     noiseMeanVar   If 'config.doReplaceWithNoise = True', you can specify the mean
                                      and variance of the Gaussian noise that will be added, by passing
                                      a tuple of (mean, variance) floats.  This overrides the
                                      "config.noiseSource" setting (but is overridden by noiseImage).
        @param[in]     references     Sequence containing reference sources detected on reference exposure.
        @param[in]     refWcs         Wcs for the reference sources.
        @param[in]     beginPriority  Only run algorithms whose 'priority' (config field that sets the order
                                      of algorithms) is greater than or equal to this value.
        @param[in]     endPriority    Only run algorithms whose 'priority' (config field that sets the order
                                      of algorithms) is less than this value.
        @return None
        """
        if references is None:
            references = [None] * len(sources)
        if len(sources) != len(references):
            raise RuntimeError("Number of sources (%d) and references (%d) don't match" %
                               (len(sources), len(references)))

        if self.config.doReplaceWithNoise and not hasattr(self, 'replaceWithNoise'):
            self.makeSubtask('replaceWithNoise')

        msg = "Measuring %d sources" % len(sources)
        if beginPriority > 0 or np.isfinite(endPriority):
            msg += ".  Measuring only priorities [%g--%g)" % (beginPriority, endPriority)
        self.log.info(msg)
        self.config.slots.setupTable(sources.table, prefix=self.config.prefix)

        # "noiseout": we will replace all the pixels within detected
        # Footprints with noise, and then add sources in one at a
        # time, measure them, then replace with noise again.  The idea
        # is that measurement algorithms might look outside the
        # Footprint, and we don't want other sources to interfere with
        # the measurements.  The faint wings of sources are still
        # there, but that's life.
        with ExposureContext(self, exposure, sources, noiseImage, noiseMeanVar), ds9.Buffering():
            # At this point the whole image should just look like noise.
            for i, (source, ref) in enumerate(zip(sources, references)):
                with SourceContext(self, exposure, sources, i):
                    # Make the measurement; note that we set refineCenter=True, but it
                    # only takes effect if config.centroider != None.
                    if ref is None:
                        self.measurer.applyWithPeak(source, exposure, True,
                                                    beginPriority, endPriority)
                    else:
                        self.measurer.applyForced(source, exposure, ref, refWcs, True,
                                                  beginPriority, endPriority)
                    if self.config.doBlendedness:
                        self.blendedness.measureChildPixels(exposure.getMaskedImage(), source)
        # Now we loop over all of the sources one more time to compute the blendedness metrics
        # on the original image (i.e. with no noise replacement).
        for source, ref in zip(sources, references):
            if self.config.doBlendedness:
                self.blendedness.measureParentPixels(exposure.getMaskedImage(), source)


    # Alias for backwards compatibility
    run = measure

class ExposureContext(object):
    """Context manager for measurement at exposure level

    At the exposure level, we replace the entire exposure with noise (if configured)
    and call the exposure-level debugging measurement hooks in the SourceMeasurementTask.
    """

    def __init__(self, task, exposure, sources, *args, **kwargs):
        """Constructor

        @param task  the SourceMeasurementTask
        @param exposure  the exposure we're measuring
        @param sources  the catalog of sources we're measuring
        @param args,kwargs  additional arguments for task.replaceWithNoise.begin
        """
        self._task = task
        self._exposure = exposure
        self._sources = sources
        self._args = args
        self._kwargs = kwargs

    def __enter__(self):
        self._task.preMeasureHook(self._exposure, self._sources)
        if self._task.config.doReplaceWithNoise:
            self._task.replaceWithNoise.begin(self._exposure, self._sources, *self._args, **self._kwargs)
        # Call the hook, with source id = -1, before we measure anything.
        # (this is *after* the sources have been replaced by noise, if noiseout)
        self._task.preSingleMeasureHook(self._exposure, self._sources, -1)

    def __exit__(self, *excInfo):
        if self._task.config.doReplaceWithNoise:
            self._task.replaceWithNoise.end(self._exposure, self._sources)
        self._task.postMeasureHook(self._exposure, self._sources)


class SourceContext(object):
    """Context manager for individual source measurement

    At the individual source level, we restore the source into the noised-out image
    (if configured) and call the source-level debugging measurement hooks in the
    SourceMeasurementTask.
    """

    def __init__(self, task, exposure, sources, index):
        """Constructor

        @param task  the SourceMeasurementTask
        @param exposure  the exposure we're measuring
        @param sources  the catalog of sources we're measuring
        @param index  the index of the source of interest in the catalog of sources
        """
        self._task = task
        self._exposure = exposure
        self._sources = sources
        self._index = index

    def __enter__(self):
        if self._task.config.doReplaceWithNoise:
            self._task.replaceWithNoise.insertSource(self._exposure, self._index)
        self._task.preSingleMeasureHook(self._exposure, self._sources, self._index)

    def __exit__(self, cls, exc, traceback):
        if self._task.config.doReplaceWithNoise:
            self._task.replaceWithNoise.removeSource(self._exposure, self._sources, self._sources[self._index])
        self._task.postSingleMeasureHook(self._exposure, self._sources, self._index)

        if exc is None:
            return

        if issubclass(cls, MemoryError):
            # always let MemoryError propagate up, as continuing just causes more problems
            return False # Don't suppress exception

        source = self._sources[self._index]
        self._task.log.warn("Error measuring source %s at %s,%s: %s" %
                            (source.getId(), source.getX(), source.getY(), exc))
        return True # Suppress exception

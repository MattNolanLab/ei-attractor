'''
Visitors related to bumps.
'''
from abc import ABCMeta, abstractmethod
import numpy as np
from scipy.optimize import leastsq
import os
import matplotlib.pyplot as plt

import analysis.spikes as aspikes
import analysis.image as image
import data_storage.sim_models.ei as simei
from interface        import DictDSVisitor 
from otherpkg.log     import getClassLogger
from analysis.image   import Position2D, fitGaussianBumpTT

from . import defaults

import logging
logger          = logging.getLogger(__name__)
speedLogger     = getClassLogger('SpeedEstimator', __name__)
speedPlotLogger = getClassLogger('SpeedPlotter', __name__)
velGainLogger   = getClassLogger('VelocityGainEstimator', __name__)


__all__ = ['BumpFittingVisitor', 'VelocityGainEstimator']

class BumpVisitor(DictDSVisitor):
    __meta__ = ABCMeta

    @abstractmethod
    def __init__(self, forceUpdate, readme, outputRoot, bumpERoot, bumpIRoot,
            winLen, tstart):
        self.forceUpdate = forceUpdate
        self.readme = readme

        self.outputRoot = outputRoot
        self.bumpERoot = bumpERoot
        self.bumpIRoot = bumpIRoot

        self.winLen = winLen
        self.tstart = tstart

    def _getSpikeTrain(self, data, monName, dimList):
        N_x = self.getNetParam(data, dimList[0])
        N_y = self.getNetParam(data, dimList[1])
        senders, times = simei.extractSpikes(data[monName])
        return senders, times, (N_x, N_y)



class BumpFittingVisitor(BumpVisitor):
    '''
    The bump fitting visitor takes spikes of the neural population (in the
    dictionary data set provided) and tries to
    fit a Gaussian shaped function to the population firing rate map (which is
    assumed to be a twisted torus). It saves the data under the 'analysis' key
    into the dataset.

    If the corresponding data is already present the visitor skips the data
    analysis and saving.
    '''
    
    def __init__(self,
            forceUpdate=False,
            readme='',
            outputRoot=defaults.analysisRoot,
            bumpERoot='bump_e',
            bumpIRoot='bump_i',
            winLen=250.0,   # ms
            tstart=None):   # ms
        super(BumpFittingVisitor, self).__init__(forceUpdate, readme,
                outputRoot, bumpERoot, bumpIRoot, winLen, tstart)

    def fitGaussianToMon(self, mon, Nx, Ny, tstart, tend):
        '''
        Fit a Gaussian function to the monitor mon, and return the results.
        '''
        senders, times = simei.extractSpikes(mon)
        sheetSize = (Nx, Ny)
        torus = aspikes.TorusPopulationSpikes(senders, times, sheetSize)
        bump = torus.avgFiringRate(tstart, tend)
        dim = Position2D(Nx, Ny)
        return fitGaussianBumpTT(bump, dim), bump


    def visitDictDataSet(self, ds, **kw):
        '''
        Apply the bump fitting procedure onto the dataset 'ds' if necessary,
        and save the data into the dataset.
        '''
        data = ds.data
        if ((self.outputRoot not in data.keys())):
            data[self.outputRoot] = {}

        a = data[self.outputRoot]

        # Determine tstart and tend
        # tstart == None --> second last time frame determined by winLen
        # tstart == full --> full simulation from the start of theta
        # tstart == int or float --> from tstart to tstart+winLen or simT
        simT = self.getOption(data, 'time')
        if (self.tstart is None):
            tstart = simT - 2*self.winLen
            tend   = tstart + winLen
        elif isinstance(self.tstart, int) or isinstance(self.tstart, float):
            tstart = self.tstart
            if winLen is not None:
                tend = tstart + self.winLen
            else:
                tend = simT
        elif (self.tstart == 'full'):
            tstart = simei.getOption(data, 'theta_start_t')
            tend = simT
        else:
            raise ValueError('Undefined tstart value: {}'.format(self.tstart))

        logger.debug("%s: tstart: %f, tend: %f", self.__class__.__name__,
                tstart, tend)
        logger.debug("%s: bumpERoot: %s", self.__class__.__name__, self.bumpERoot)

        if (self.bumpERoot not in a.keys()) or self.forceUpdate:
            logger.info("%s: Analysing an E dataset", self.__class__.__name__)
            # Fit the Gaussian onto E neurons
            Nx  = self.getNetParam(data, 'Ne_x')
            Ny  = self.getNetParam(data, 'Ne_y')
            mon = data['spikeMon_e']
            ((A, mu_x, mu_y, sigma), err2), bump = self.fitGaussianToMon(mon,
                    Nx, Ny, tstart, tend)
            a[self.bumpERoot] = {
                    'A' : A,
                    'mu_x' : mu_x,
                    'mu_y' : mu_y,
                    'sigma' : np.abs(sigma),
                    'err2'  : np.sum(err2),
                    'bump_e_rateMap' : bump
            }
        else:
            logger.info("%s: E data present. Skipping analysis.",
                    self.__class__.__name__)

        # Only export population firing rates of I neurons
        if (self.bumpIRoot not in a.keys() or self.forceUpdate):
            logger.info("%s: Analysing an I dataset", self.__class__.__name__)
            Nx  = simei.getNetParam(data, 'Ni_x')
            Ny  = simei.getNetParam(data, 'Ni_y')
            senders, times = simei.extractSpikes(data['spikeMon_i'])
            sheetSize = (Nx, Ny)
            torus = aspikes.TorusPopulationSpikes(senders, times, sheetSize)
            bump = torus.avgFiringRate(tstart, tend)
            a[self.bumpIRoot] = dict(
                    bump_i_rateMap=bump)
        else:
            logger.info("%s: I data present. Skipping analysis.",
                    self.__class__.__name__)


###############################################################################


def fitCircularSlope(bumpPos, times, normFac):
    '''
    Fit a (circular) slope line onto velocity response of the bump and extract
    the slope (velocity), in neurons/s. The start position of the bump is
    ignored here (i.e. bumpPos -= bumpPos[0])

    Parameters
    ----------
    bumpPos : numpy array
        An array of bump positions
    times : numpy array
        A corresponding vector of times, of the same size as bumpPos
    normFac : float
        Normalizing factor for the positional vector. In fact this is the size
        of the toroidal sheet (X direction).
    output : float
        Estimated bump speed (neurons/time unit).
    '''
    bumpPos = np.asarray(bumpPos)
    t = np.array(times) - times[0]
    bumpPos_norm = np.unwrap(1.0 * bumpPos / normFac * 2 * np.pi) # normalise to 2*Pi
    bumpPos_norm -= bumpPos_norm[0]
    func = lambda X: X[0]*t - bumpPos_norm
    x0 = np.array([0.0])  # slope
    x = leastsq(func, x0)
    return x[0][0] / 2. / np.pi * normFac
    

def getLineFit(X, Y):
    '''
    Fit a line to data
    '''
    func = lambda P: P[0]*X  - Y
    P0 = np.array([0.0]) # slope
    P = leastsq(func, P0)
    return P[0][0]*X, P[0][0]


def fitBumpSpeed(IvelVec, bumpSpeed, fitRange):
    '''
    Fit a line to bumpSpeed vs IvelVec, useing IvelVec[0:fitRange+1]
    '''
    fitSpeed    = np.array(bumpSpeed)[:, 0:fitRange+1]
    fitIvelVec  = np.repeat([IvelVec[0:fitRange+1]], fitSpeed.shape[0],
            axis=0)
    fitSpeed    = fitSpeed.flatten()
    fitIvelVec   = fitIvelVec.flatten()

    # Get rid of NaNs here
    nanIdx = np.isnan(fitSpeed)
    fitSpeed = fitSpeed[np.logical_not(nanIdx)]
    fitIvelVec = fitIvelVec[np.logical_not(nanIdx)]

    line, slope = getLineFit(fitIvelVec, fitSpeed)
    lineFitErr  = np.abs(line - fitSpeed) / len(fitIvelVec)

    # Compose return values; as a function of IvelVec[0:fitRange+1]
    retLine = np.array(IvelVec[0:fitRange+1]) * slope
    retIvelVec = IvelVec[0:fitRange+1]
    return retLine, slope, lineFitErr, retIvelVec


def nanPaddedArray(data, size=None):
    r = size[0] if size is not None else None
    c = size[1] if size is not None else None
    if r is None:
        r = len(data)
    if c is None:
        c = 0
        for row in data:
            c = len(row) if len(row) > c else c

    res = np.ndarray((r, c), dtype=np.double) * np.nan
    for rowIdx, row in enumerate(data):
        res[rowIdx, 0:len(row)] = row

    return res


class SpeedEstimator(BumpVisitor):
    '''Estimate raw relationship between injected (velocity) current and bump
    speed.

    Steps performed:

        1. Do the estimation.

        2. Store the raw data for each velocity current value.

        3. Store the data into a top-level stracture in the data hierarchy
    '''
    allowedAxis = ['horizontal', 'vertical']

    def __init__(self,
            win_dt=20.0,
            winLen=250.0,
            forceUpdate=False,
            outputRoot=defaults.analysisRoot,
            axis='vertical',
            changeSign=False,
            readme=''):
        super(SpeedEstimator, self).__init__(
                forceUpdate,
                readme,
                outputRoot,
                None,
                None,
                winLen,
                None)
        self.win_dt = win_dt
        self.axis = self.checkAxis(axis)
        self.changeSign = changeSign

    def checkAxis(self, axis):
        if axis not in self.allowedAxis:
            raise ValueError("axis parameter must be one of {0}".format(self.allowedAxis))
        return axis

    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        trials = data['trials']

        slopes = []
        for trialNum in xrange(len(trials)):
            speedLogger.info("Trial no. %d/%d", trialNum, len(trials)-1)
            slopes.append([])
            IvelVec = trials[trialNum]['IvelVec']
            for IvelIdx in xrange(len(IvelVec)):
                iData = trials[trialNum]['IvelData'][IvelIdx]
                speedLogger.info("Processing vel. index: %d/%d; %.2f pA",
                        IvelIdx, len(IvelVec) - 1, IvelVec[IvelIdx])

                slope = None
                if self.outputRoot in iData.keys() and not self.forceUpdate:
                    speedLogger.info("Data present. Adding slope to the list.")
                    slope = iData[self.outputRoot]['slope']
                else:
                    speedLogger.debug("\tEstimating bump positions...")
                    senders, times, sheetSize =  self._getSpikeTrain(iData,
                            'spikeMon_e', ['Ne_x', 'Ne_y'])
                    pop = image.SingleBumpPopulation(senders, times, sheetSize)
                    tStart = self.getOption(iData, 'theta_start_t')
                    tEnd   = self.getOption(iData, 'time') - 2*self.win_dt
                    fitList = pop.bumpPosition(tStart, tEnd, self.win_dt,
                            self.winLen, fullErr=False)

                    if self.axis == 'vertical':
                        bumpPos = fitList.mu_y
                        axisSize = pop.Ny
                    elif self.axis == 'horizontal':
                        bumpPos = fitList.mu_x
                        axisSize = pop.Nx

                    speedLogger.debug("\tFitting the circular slope...")
                    # NOTE: If the bump moves in an opposite direction, we have to
                    # negate the speed
                    slope = fitCircularSlope(bumpPos, fitList.times, axisSize)
                    slope *= 1e3 # Correction msec --> sec
                    if self.changeSign:
                        slope *= -1
                    speedLogger.debug("\tslope: {0:.3f}".format(slope))

                    speedLogger.debug("\tSaving data for the current velocity index.")
                    iData[self.outputRoot] = dict(
                            positions = dict(
                                A            = np.asarray(fitList.A),
                                mu_x         = np.asarray(fitList.mu_x),
                                mu_y         = np.asarray(fitList.mu_y),
                                sigma        = np.asarray(fitList.sigma),
                                err2         = np.asarray(fitList.err2),
                                ln_L         = np.asarray(fitList.ln_L),
                                lh_precision = np.asarray(fitList.lh_precision),
                                times        = np.asarray(fitList.times)
                            ),
                            slope = slope
                    )
                    iData.flush()

                assert(slope is not None)
                slopes[trialNum].append(slope)

        slopes = nanPaddedArray(slopes)
        analysisTop = {'bumpVelAll' : slopes}
        data[self.outputRoot] = analysisTop

        printoptions_orig = np.get_printoptions()
        np.set_printoptions(precision=3, threshold=np.infty, linewidth=100,
                            suppress=True)
        speedLogger.info("Slopes:\n%s", slopes)
        np.set_printoptions(**printoptions_orig)


def plotVelocity(IvelVec, slopes, ax):
    # Plot the estimated bump velocities (nrns/s)
    avgSlope    = np.nanmean(slopes, axis = 0)
    nTrials     = slopes.shape[0]
    stdErrSlope = np.nanstd(slopes, axis  = 0) / np.sqrt(nTrials)

    ax.errorbar(IvelVec, avgSlope, stdErrSlope, fmt='o-')
    ax.plot(IvelVec, slopes.T, 'o', color='blue', alpha=0.4)
    ax.set_xlabel('Velocity current (pA)')
    ax.set_ylabel('Bump velocity (neurons/s)')

    leg = ['Average bump speed', 'Estimated bump speed']

    return leg


class SpeedPlotter(DictDSVisitor):
    def __init__(self,
                 rootDir=None,
                 plotFittedLine=True,
                 readme=''):
        '''
        Parameters
        ----------
        rootDir : str
            Root output directory where to save the image. It will be a
            subdirectory of where the data is located.
        '''
        self.rootDir = rootDir
        self.plotFittedLine = plotFittedLine
        self.readme = readme


    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        trials = data['trials']
        IvelVec = trials[0]['IvelVec'] # All the same
        a      = data[defaults.analysisRoot]
        slopes = a['bumpVelAll']
        

        if 'fileName' not in kw.keys():
            msg = 'Did not receive the fileName as a keyword argument.'
            velGainLogger.warn(msg)
            return

        speedPlotLogger.info("Plotting the velocity points")
        fig = plt.figure()
        ax = fig.gca()
        leg = plotVelocity(IvelVec, slopes, ax)

        if self.plotFittedLine:
            speedPlotLogger.info("Plotting the fitted line")
            ax.plot(a['fitIvelVec'], a['lineFitLine'], 'o-')
            t = "Line fit slope: {0:.3f} nrns/s/pA, error: {1:.3f} " + \
                "neurons/s (norm)"
            ax.set_title(t.format(a['lineFitSlope'], np.sum(a['lineFitErr'])))
            leg.append('Line fit')

        ax.legend(leg, loc='best')
        
        dataRootDir, dataFileName = os.path.split(kw['fileName'])
        figDir = dataRootDir + "/"
        if self.rootDir is not None:
            figDir += self.rootDir
            if not os.path.exists(figDir):
                os.makedirs(figDir)

        figFileName = os.path.splitext(dataFileName)[0] + '.pdf'
        figPath = "{0}/{1}".format(figDir, figFileName)
        speedPlotLogger.info("Saving figure to '{0}'".format(figPath))
        fig.savefig(figPath)


class VelocityGainEstimator(BumpVisitor):
    '''
    A visitor that estimates the relationship between injected velocity current
    and bump speed.
    '''
    def __init__(self,
            bumpSpeedMax,
            forceUpdate=False,
            outputRoot=defaults.analysisRoot,
            readme=''):
        super(VelocityGainEstimator, self).__init__(
                forceUpdate,
                readme,
                outputRoot,
                None,
                None,
                None,
                None)
        assert(bumpSpeedMax is not None)
        self.bumpSpeedMax = bumpSpeedMax # cm/s


    def visitDictDataSet(self, ds, **kw):
        '''
        Visit a data set that contains all the trials and Ivel simulations.
        '''
        data = ds.data
        trials = data['trials']
        slopes = data[defaults.analysisRoot]['bumpVelAll']
        IvelVec = trials[0]['IvelVec'] # All the same

        velGainLogger.info('Fitting lines')
        if (len(trials) == 0):
            velGainLogger.warn('Something wrong: len(trials) == 0. Not plotting data')
            return

        # Fit a line (nrns/s/pA)
        # results
        line       = None
        lineFitErr = None
        slope      = None
        fitIvelVec = None
        errSum     = None
        # All, in case we need the one with max range
        line_all       = []
        lineFitErr_all = []
        slope_all      = []
        fitIvelVec_all = []
        maxRange_all   = []
        for fitRange in xrange(1, len(IvelVec)):
            velGainLogger.info("fitRange: {0}".format(fitRange))
            newLine, newSlope, newErr, newIvelVecRange = \
                    fitBumpSpeed(IvelVec, slopes, fitRange)
            line_all.append(newLine)
            lineFitErr_all.append(newErr)
            slope_all.append(newSlope)
            fitIvelVec_all.append(newIvelVecRange)
            maxRange_all.append(newLine[-1])
            # Keep only lines that covers the desired range and have a
            # minimal error of fit
            if ((newLine[-1] >= self.bumpSpeedMax)):
                errSumNew = np.sum(newErr)
                if ((line is None) or (errSumNew <= errSum)):
                    line       = newLine
                    lineFitErr = newErr
                    slope      = newSlope
                    fitIvelVec = newIvelVecRange
                    errSum     = errSumNew

        if (line is None):
            velGainLogger.info(\
                ('No suitable fits that cover <0, {0:.2f}> neurons/s.' +\
                ' Using the fit with the max. bump speed range.').format(\
                self.bumpSpeedMax))
            velGainLogger.info("Bump speed maxima: {0}".format(maxRange_all))
            maxRangeIdx = np.argmax(maxRange_all)
            line        = line_all[maxRangeIdx]
            lineFitErr  = lineFitErr_all[maxRangeIdx]
            slope       = slope_all[maxRangeIdx]
            fitIvelVec  = fitIvelVec_all[maxRangeIdx]

        velGainLogger.info('Saving the line fit data (top level)')
        velGainLogger.info(\
                "Estimated bump velocity gain: %.3f nrns/s/pA, err=%.3f nrns/s",
                slope, np.sum(lineFitErr))

        data[self.outputRoot].update({
            'lineFitLine'  : line,
            'lineFitSlope' : slope,
            'lineFitErr'   : lineFitErr,
            'fitIvelVec'   : fitIvelVec
        })



##############################################################################
# Bump position visitor
class BumpPositionVisitor(BumpVisitor):
    '''
    A visitor that estimates position of the bump attractor by fitting a
    Gaussian function onto the population vector.

    The population firing is sliced onto windows of len winLen
    '''
    def __init__(self,
            forceUpdate=False,
            readme='',
            outputRoot=defaults.analysisRoot,
            bumpERoot='bump_e',
            bumpIRoot='bump_i',
            tstart=None,        # ms
            tend=None,          # ms
            winLen=250.0,       # ms
            win_dt=25.0,        # ms
            rateThreshold=.1,   # Hz
            maxRelBumpDiam=1.0): # fraction of the torus size
        super(BumpPositionVisitor, self).__init__(forceUpdate, readme,
                outputRoot, bumpERoot, bumpIRoot, winLen, tstart)
        self.tend = tend
        self.win_dt = win_dt
        self.rateThreshold  = rateThreshold
        self.maxRelBumpDiam = maxRelBumpDiam
        

    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        data.setItemChained((self.outputRoot, self.bumpERoot), {},
                overwriteLast=False)
        out = data[self.outputRoot][self.bumpERoot]

        # Resolve times
        if self.tstart is None:
            tstart = self.getOption(data, 'theta_start_t')
        else:
            tstart = self.tstart
        if self.tend is None:
            tend = self.getOption(data, 'time') - self.win_dt

        senders, times, sheetSize =  self._getSpikeTrain(data, 'spikeMon_e',
                ['Ne_x', 'Ne_y'])
        pop = image.SingleBumpPopulation(senders, times, sheetSize)

        # Bump fits
        if 'positions' not in out.keys() or self.forceUpdate:
            logger.info('%s: Analysing data set (pos).', self.__class__.__name__)
            bumpFits = pop.bumpPosition(tstart, tend, self.win_dt, self.winLen,
                    fullErr=False)
            out['positions'] = dict(
                    A            = np.asarray(bumpFits.A),
                    mu_x         = np.asarray(bumpFits.mu_x),
                    mu_y         = np.asarray(bumpFits.mu_y),
                    sigma        = np.asarray(bumpFits.sigma),
                    err2Sum      = np.asarray(bumpFits.err2),
                    ln_L         = np.asarray(bumpFits.ln_L),
                    lh_precision = np.asarray(bumpFits.lh_precision),
                    times        = np.asarray(bumpFits.times),
                    readme       = self.readme
            )
        else:
            logger.info('%s: Positional data already present. Skipping analysis.',
                    self.__class__.__name__)

        # Uniform fits (ML)
        if 'uniformML' not in out.keys() or self.forceUpdate:
            logger.info('%s: Analysing data set (uni).', self.__class__.__name__)
            uniformFits = pop.uniformFit(tstart, tend, self.win_dt,
                    self.winLen, fullErr=False)
            out['uniformML'] = dict(
                    mu     = np.asarray(uniformFits.mu),
                    sigma2 = np.asarray(uniformFits.sigma2),
                    ln_L   = np.asarray(uniformFits.ln_L),
                    err2   = np.asarray(uniformFits.err2),
                    times  = np.asarray(uniformFits.times)
            )
        else:
            logger.info('%s: Uniform data already present. Skipping analysis.',
                    self.__class__.__name__)

        
        sigma = out['positions/sigma']
        A     = out['positions/A']
        bumpDiam = sigma * np.sqrt(-2. * np.log(self.rateThreshold / A))
        isBumpFrames = bumpDiam < self.maxRelBumpDiam * np.min(sheetSize)
        nFrames = len(isBumpFrames)
        out['isBump'] = dict(
                isBumpFrames = isBumpFrames,
                fracTotal = np.count_nonzero(isBumpFrames) / float(nFrames)
        )




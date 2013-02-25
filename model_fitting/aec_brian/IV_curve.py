from brian.library.electrophysiology import *
import numpy as np
from scipy.optimize import fmin, leastsq


pF = 1e-12

def voltage_bins(V, binStart, binEnd, nbins):
    '''
    Picks indices of V that are within bin ranges.
    Samples outside binStart and binEnd are discarded
    '''
    if binStart > binEnd:
        raise Exception("binStart must be less than binEnd")
    binCenters = np.linspace(binStart, binEnd, nbins+1)
    binWidth = (binEnd - binStart)/nbins
    Ibin = []
    for c in binCenters:
        Ibin.append((np.logical_and(V > c-binWidth/2, V <=
            c+binWidth/2)).nonzero()[0]) 

    return (binCenters, Ibin)
    

def findIcurveZeroCrossings(V, I):
    '''
    Take scatter data of V and I and find where they have left and right zero
    crossings
    '''
    binStart = -80e-3
    binEnd   = -35e-3
    nbins    = 25
    binCenters, bin_it = voltage_bins(V, binStart, binEnd, nbins)
    # Find first and second zero crossings. If the curve is reasonably nice,
    # we will find only two of them

    

class Kernels(object):
    ''' Structure to hold kernels'''
    def __init__(self, Kfull=None, Ke=None, Km=None):
        self.Kfull = Kfull
        self.Ke = Ke
        self.Km = Km


def compensate_voltage(Iin, Vrec, Iin_ker, Vm_ker, ksize=15e-3, tail_start=3e-3):
    '''
    Estimate membrane and electrode kernel and substract the electrode kernel
    filtered voltage from Vm.
    '''
    Kfull, Vmean = full_kernel(Vm_ker, Iin_ker, ksize, full_output=True)
    Ke, Km = electrode_kernel(Kfull, tail_start, full_output=True)
    Vcomp = AEC_compensate(Vrec, Iin, Ke)
    return Vcomp, Kernels(Kfull, Ke, Km), Vmean


class DynamicIVCurve(object):
    '''
    A basic class that generically estimates an IV curve from a given set of
    recorded membrane voltage, input current and the derivative of the same
    recorded membrane voltage (multiplied by dt).

    It uses Brette et al. Active Electrode Compensation module implemented in
    brian (done offline). Therefore brian module is required in fact for this.

    This class doesn't do any assumptions about what time indices it should use.
    '''

    class Curve:
        def __init__(self, V, I, Istd, IN, C):
            '''
            V       V binned
            I       I binned
            Istd    Std. deviation of Ibins
            IN      Number of points in each Ibins
            C       Estimated membrane capacitance
            '''
            self.V      = V
            self.I      = I
            self.Istd   = Istd
            self.IN     = IN
            self.C      = C

    class ExpIaFParams:
        '''
        Parameters of an exp. IaF model fit
        '''
        def __init__(self, C, taum, Em, dT, VT):
            self.C = C
            self.taum = taum
            self.Em = Em
            self.dT = dT
            self.VT = VT

        def getCurve(self, V):
            return 1/self.taum*(self.Em - V + self.dT*np.exp((V - self.VT)/self.dT))

        def __repr__(self):
            ret  = "C:    " + str(self.C*1e12) + " pF\n"
            ret += "taum: " + str(self.taum*1e3) + " ms\n"
            ret += "Em:   " + str(self.Em*1e3) + " mV\n"
            ret += "dT:   " + str(self.dT*1e3) + " mV\n"
            ret += "VT:   " + str(self.VT*1e3) + " mV\n"
            return ret

    #def binCurrent(self, Vm, Im, 

    def __init__(self, Iin, Vm, dV, dt, binStartV, binEndV, nbins, C=None):
        self.Iin = Iin
        self.Vm = Vm
        self.dV = dV
        self.binStartV = binStartV
        self.binEndV = binEndV
        self.nbins = nbins
        self.dt = dt
        self._expFit = None

        # Estimate voltage bin range

        self.binCenters, self.binIds = voltage_bins(Vm, binStartV, binEndV, nbins)
        
        self.Cest_bin_id = 10
        if C is None:
            self.C0 = 300*pF
            self.Cest = self.findCapacitance(Iin[self.binIds[self.Cest_bin_id]],
                    dV[self.binIds[self.Cest_bin_id]]/dt, self.C0)
            print("Estimated capacitance: " + str(self.Cest/pF) + " pF.")
        else:
            self.Cest = C
        self.Im = Iin - self.Cest*dV/dt

        self.Imean = np.zeros(len(self.binCenters))
        self.Istd  = np.zeros(len(self.binCenters))
        self.IN    = np.zeros(len(self.binCenters))
        for i in xrange(len(self.binCenters)):
            self.Imean[i] = np.mean(self.Im[self.binIds[i]])
            self.Istd[i]  = np.std(self.Im[self.binIds[i]])
            self.IN[i]    = len(self.binIds[i])

        self.cur = DynamicIVCurve.Curve(self.binCenters, self.Imean, self.Istd, self.IN, self.Cest)


    def findCapacitance(self, Iin, dVdt, C0):
        '''
        Estimate membrane capacitance, by using input current and dV/dt at a
        specific voltage value (set by the user, ideally somewhere near resting Vm)
        '''
        return fmin(lambda Ce: np.var(Iin/Ce - dVdt), C0)[0]

    def fitExponentialIaF(self):
        '''Fit an exponential integrate and fire model to the data'''
        if self._expFit is None:
            taum0 = 10e-3
            Em0 = -80e-3
            dT0 = 2e-3
            VT0 = -55e-3
            x0 = np.array([taum0, Em0, dT0, VT0])
            FV = -self.Im/self.cur.C

            Vrange_it   = np.logical_and(self.Vm >= self.binCenters[0], self.Vm <= self.binCenters[-1])
            Vrange      = self.Vm[Vrange_it]
            FV_range    = FV[Vrange_it]

            fun = lambda x: 1/x[0]*(x[1] - Vrange + x[2]*np.exp((Vrange -
                x[3])/x[2])) - FV_range
            xest, ierr  = leastsq(fun, x0, full_output=0)
            self._expFit = DynamicIVCurve.ExpIaFParams(C=self.cur.C, taum=xest[0], Em=xest[1], dT=xest[2], VT=xest[3])
            return self._expFit
        else:
            return self._expFit



class DynamicIVCurveAfter(DynamicIVCurve):
    '''
    This class computes the dynamic IV curve, but excludes all the point which
    are < tafter seconds after the peak of previous spikes.

    Input currents and voltages must be contiguous regions.
    It is assumed that sampling rate is constant (times[i+1] - times[i] = dt)
    '''
    def __init__(self, Iin, Vm, times, binStartV, binEndV, nbins, tafter,
            C=None, Vth=0):

        self.Iin_all = Iin
        self.Vm_all  = Vm
        self.times_all = times
        self.tafter = tafter
        self.Vth = Vth
        self.dt = times[1] - times[0]

        time_id = self.pickVPreSpike(Vm, tafter, Vth)
        time_id = time_id[0:len(time_id)-1]  # If last spike is the last in the array

        DynamicIVCurve.__init__(self, Iin[time_id], Vm[time_id], np.diff(Vm)[time_id],
                self.dt, binStartV, binEndV, nbins, C)
        self.time_id = time_id
        self.times = self.times_all[time_id]


    def pickVPreSpike(self, Vm, t_after, V_th):
        '''
        Pick only voltage traces that are more than t_after time units after the
        last spike. By setting V_th, one can fine_tune spike threshold.
        '''
        spike_times = spike_peaks(Vm, V_th)
        time_ids = np.arange(len(Vm))
        t_mask = np.ndarray(len(Vm), dtype=bool)
        t_mask[:] = True
        for t_spike in spike_times:
            t_mask = np.logical_and(t_mask, np.logical_or(time_ids <= t_spike,
                time_ids > t_spike+t_after))
    
        return time_ids[t_mask]
    
class DynamicIVCurveAfterRegion(DynamicIVCurve):
    '''
    Compute the dynamic IV curve from bounded regions after the spike.
    '''
    def __init__(self, Iin, Vm, times, binStartV, binEndV, nbins, region,
            C=None, Vth=0):
        self.Iin_all = Iin
        self.Vm_all = Vm
        self.times_all = times
        self.Vth = Vth
        self.dt = times[1] - times[0]
        self.region = region
        self._region = (int(region[0]/self.dt), int(region[1]/self.dt))
        if region[0] >= region[1]:
            raise Exception('region[0] < region[1] !')

        self.time_id = self.pickVAfterSpikeRegion(Vm, self._region, Vth)
        self.time_id = self.time_id[0:len(self.time_id)-1]

        DynamicIVCurve.__init__(self, Iin[self.time_id], Vm[self.time_id],
                np.diff(Vm)[self.time_id], self.dt, binStartV, binEndV, nbins, C)
        self.times = self.times_all[self.time_id]

    def pickVAfterSpikeRegion(self, Vm, region, Vth):
        '''Pick only a regions of length region[1] - region[0] after the peak of
        each spike in the trace'''
        spike_times = spike_peaks(Vm, Vth)
        time_ids = np.arange(len(Vm))
        t_mask = np.ndarray(len(Vm), dtype=bool)
        t_mask[:] = False
        for t_spike in spike_times:
            t_mask = np.logical_or(t_mask, np.logical_and(time_ids >=
                t_spike+region[0], time_ids < t_spike+region[1]))
            t_mask = np.logical_and(t_mask, np.logical_or(time_ids <= t_spike,
                time_ids >= t_spike+region[0]))

        return time_ids[t_mask]


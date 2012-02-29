from brian.library.electrophysiology import *
from numpy import diff
from scipy.optimize import fmin


pF = 1e-12

class Kernels(object):
    ''' Structure to hold kernels'''
    def __init__(self, Kfull=None, Ke=None, Km=None):
        self.Kfull = Kfull
        self.Ke = Ke
        self.Km = Km


def compensate_voltage(Iin, Vrec, Vm_ker, Iin_ker, ksize=15e-3, tail_start=3e-3):
    '''
    Estimate membrane and electrode kernel and substract the electrode kernel
    filtered voltage from Vm.
    '''
    K_full, Vmean = full_kernel(Vm_ker, Im_ker, ksize, full_output=True)
    Ke, Km = electrode_kernel(K_full, tail_start, full_output=True)
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

    def __init__(self, Iin, Vm, dV, dt, binStartV, binEndV, nbins):
        self.Iin = Iin
        self.Vm = Vm
        self.dV = dV
        self.binStartV = binStartV
        self.binEndV = binEndV
        self.nbins = nbins
        self.dt = dt

        self.binCenters, self.binIds = voltage_bins(Vm, binStartV, binEndV, nbins)
        
        self.Cest_bin_id = 27
        self.C0 = 300*pF
        self.Cest = findCapacitance(Iin[self.binIds[self.Cest_bin_id]],
                dV[self.binIds[self.Cest_bin_id]]/dt, self.C0)
        print("Estimated capacitance: " + str(self.Cest/pF) + " pF.")
        self.Im = self.Cest*dV/dt - Iin

        self.Imean = np.zeros(len(binCenters))
        self.Istd  = np.zeros(len(binCenters))
        self.IN    = np.zeros(len(binCenters))
        for i in xrange(len(binCenters)):
            self.Imean[i] = np.mean(Im[self.binIds[i]])
            self.Istd[i]  = np.std(Im[self.binIds[i]])
            self.IN[i]    = len(self.binIds[i])


    def voltage_bins(self, V, binStart, binEnd, nbins):
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
    
    
    def findCapacitance(self, Iin, dVdt, C0):
        '''
        Estimate membrane capacitance, by using input current and dV/dt at a
        specific voltage value (set by the user, ideally somewhere near resting Vm)
        '''
        return fmin(lambda Ce: np.var(Iin/Ce - dVdt), C0, xtol=0.00001)[0]


class DynamicIVCurveAfter(DynamicIVCurve):
    '''
    This class computes the dynamic IV curve, but excludes all the point which
    are < tafter seconds after the peak of previous spikes.

    Input currents and voltages must be contiguous regions.
    It is assumed that sampling rate is constant (times[i+1] - times[i] = dt)
    '''
    def __init__(self, Iin, Vm, times, binStartV, binEndV, nbins, tafter, Vth=0):

        self.Iin_all = Iin
        self.Vm_all  = Vm
        self.times = times
        self.tafter = tafter
        self.Vth = Vth
        self.dt = times[1] - times[0]

        time_id = self.pickVPreSpike(Vm, tafter)
        time_id = time_id[0:len(time_id)-1]  # If last spike is the last in the array

        DynamicIVCurve.__init__(Iin[time_id], Vm[time_id], np.diff(Vm)[time_id],
                dt, binStartV, binEndV, nbins)


    def pickVPreSpike(self, Vm, t_after, V_th=0):
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
            Vth=0):
        pass
    


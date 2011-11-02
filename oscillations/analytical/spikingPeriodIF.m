function [T I] = spikingPeriodIF(opt, refrac, I0, Imax)
    
    nPoints = 100;

    I = linspace(I0, Imax, nPoints);
    T = opt.taum * log(opt.Rm * I ./ (opt.Rm*I - opt.Vt + opt.Vr)) + refrac;

end
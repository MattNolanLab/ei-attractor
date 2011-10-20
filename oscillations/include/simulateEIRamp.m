function [spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEIRamp(o, net_data)
    % Simulation of excitatory-inhibitory neural loop
    % Spikerecord_e/i -- sparse (every timestep) spiking activity
    % spikeTimes - dense activity sorted out by spike time
    
    
    % All variables are in basic units, i.e. s, volt, etc.
    Ne = o.Ne;
    Ni = o.Ni;

    % Excitatory cells
    taum_e = o.taum_e;
    taue = o.taue;
    El_e = o.El_e;
    Vt_e = o.Vt_e;
    Vr_e = o.Vr_e;
    e_sparseness = o.e_sparseness;
    Ie = o.Ie;
    we = o.we;


    % Inhibitory cell
    taum_i = o.taum_i;
    taui = o.taui;
    El_i = o.El_i;
    Vt_i = o.Vt_i;
    Vr_i = o.Vr_i;
    i_sparseness = o.i_sparseness;
    Ii = o.Ii;
    wi = o.wi;

    % Noise normalized per time unit (ms)
    noise_sigma = o.noise_sigma;

    % Euler settings
    dt = o.dt;
    % Times
    T = o.T;
    times = 0:dt:T;


    % Build excitatory neurons state
%     Ve = Vr_e + (Vt_e - Vr_e) * rand(Ne, 1);
%     Vi = Vr_i + (Vt_i - Vr_i) * rand(Ni, 1);

    Ve = Vr_e + o.sigma_init_cond * rand(Ne, 1);
    Vi = Vr_i + o.sigma_init_cond * rand(Ni, 1);


    ge = zeros(Ni, 1);
    gi = zeros(Ne, 1);

    spikeMon_e = {};
    spikeMon_i = {};

    Emon_i = o.Emon_i;
    Imon_i = o.Imon_i;
    Vmon_e = zeros(numel(Emon_i), numel(times));
    Vmon_i = zeros(numel(Imon_i), numel(times));
    Vmon_ge = zeros(numel(Imon_i), numel(times));
    Vmon_gi = zeros(numel(Emon_i), numel(times));
    Vmon_t = times;
    
    % Setup connections Mij: j --> i
    % Assuming constant and uniform excitatory/inhibitory weights
    Me = net_data.Me;
    Mi = net_data.Mi;

    spikeRecord_e = sparse(Ne, size(times, 2));
    spikeRecord_i = sparse(Ni, size(times, 2));

%    display 'Simulation running...'
    t = 0;

    fired_e = zeros(Ne, 1);
    fired_i = zeros(Ni, 1);
    f_i = 1;
    t_spike = [];

    t_i = 1;
    for t = times
        fired_e = Ve > Vt_e;
        fired_i = Vi > Vt_i;

        Ve(fired_e) = Vr_e;
        Vi(fired_i) = Vr_i;

        spikeRecord_e(:, t_i) = double(fired_e);
        spikeRecord_i(:, t_i) = double(fired_i);
        
        gi = gi + Mi*fired_i * wi;
        ge = ge + Me*fired_e * we;

        
        Vmon_e(:, t_i) = Ve(Emon_i);
        Vmon_i(:, t_i) = Vi(Imon_i);
        Vmon_e(fired_e(Emon_i), t_i) = o.spikeVm;
        Vmon_i(fired_i(Imon_i), t_i) = o.spikeVm;
        Vmon_ge(:, t_i) = ge(Imon_i);
        Vmon_gi(:, t_i) = gi(Emon_i);


        % Check if neurons fired and add to syn. conductances
        dVe = dt * 1/taum_e * (El_e - Ve - gi + Ie);
        dVi = dt * 1/taum_i * (El_i - Vi + ge + Ii);

        dge = dt * -1/taue * ge;
        dgi = dt * -1/taui * gi;

        Ve = Ve + dVe + noise_sigma*randn(Ne, 1);
        Vi = Vi + dVi + noise_sigma*randn(Ni, 1);
        ge = ge + dge;
        gi = gi + dgi;
        
        Ie = Ie + o.dIe*dt;
        Ii = Ii + o.dIi*dt;

        t_i = t_i + 1;
    end


    Vmon.e = Vmon_e;
    Vmon.i = Vmon_i;
    Vmon.ge = Vmon_ge;
    Vmon.gi = Vmon_gi;
    Vmon.t = Vmon_t;

end
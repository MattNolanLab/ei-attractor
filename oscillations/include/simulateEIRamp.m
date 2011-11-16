function [spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEIRamp(o, net_data)
    % Simulation of excitatory-inhibitory neural loop
    % Spikerecord_e/i -- sparse (every timestep) spiking activity
    % spikeTimes - dense activity sorted out by spike time
    
    
    % All variables are in basic units, i.e. s, volt, etc.
    Ne = o.Ne;
    Ni = o.Ni;
    
    % Euler settings
    dt = o.dt;
    % Times
    T = o.T;
    times = 0:dt:T;



    % Excitatory cells
    taum_e = o.taum_e;
    taue = o.taue;
    El_e = o.El_e;
    Vt_e = o.Vt_e;
    Vr_e = o.Vr_e;
    Rm_e = o.Rm_e;
    e_sparseness = o.e_sparseness;
    Ie = o.Ie;
    we = o.we;
    tau_adapt_e = o.refrac_e;


    % Inhibitory cell
    taum_i = o.taum_i;
    taui_dec = o.taui_dec;
    taui_rise = o.taui_rise;
    El_i = o.El_i;
    Vt_i = o.Vt_i;
    Vr_i = o.Vr_i;
    Rm_i = o.Rm_i;
    i_sparseness = o.i_sparseness;
    Ii = o.Ii;
    wi = o.wi;
    tau_adapt_i = o.refrac_i;
    
    V_rev_i = o.V_rev_i;
    V_rev_e = o.V_rev_e;
    
    Vclamp = o.Vclamp;

    % Noise normalized per time unit (ms)
    noise_sigma = o.noise_sigma;


    % Build excitatory neurons state
%     Ve = Vr_e + (Vt_e - Vr_e) * rand(Ne, 1);
%     Vi = Vr_i + (Vt_i - Vr_i) * rand(Ni, 1);

    Ve = Vr_e + o.sigma_init_cond * rand(Ne, 1);
    Vi = Vr_i + o.sigma_init_cond * rand(Ni, 1);


    ge = zeros(Ni, 1);
    gi1 = zeros(Ne, 1);
    gi2 = zeros(Ne, 1);
    
    % Adaptation constants
    g_adapt_e = zeros(Ne, 1);
    g_adapt_i = zeros(Ni, 1);

    spikeMon_e = {};
    spikeMon_i = {};

    Emon_i = o.Emon_i;
    Imon_i = o.Imon_i;
    Vmon_e = zeros(numel(Emon_i), numel(times));
    Vmon_i = zeros(numel(Imon_i), numel(times));
    Vmon_Isyn_e = zeros(numel(Emon_i), numel(times));
    Vmon_Isyn_i = zeros(numel(Imon_i), numel(times));
    Vmon_t = times;
    
    % Setup connections Mij: j --> i
    % Assuming constant and uniform excitatory/inhibitory weights
    Me = net_data.Me;
    Mi = net_data.Mi;
    
    % Draw synaptic weights from exponential distributions with mean given
    % by we/wi
    
    % Excitatory connections
    for we_it = 1:Ni
        %Me(we_it, :) = Me(we_it, :) .* exprnd(o.we, Ne);
        e_V = o.we_std^2;
        e_MU = log(o.we^2 / sqrt(e_V+o.we^2));
        e_SIGMA = sqrt(log(e_V/o.we^2 + 1));

        Me(we_it, :) = Me(we_it, :) .* lognrnd_local(e_MU, e_SIGMA, 1, Ne);
    end
%     % Inhibitory_weights
%     for wi_it = 1:Ne
%         Mi(wi_it, :) = Mi(wi_it, :) .* sampleExponential(o.wi, Ni);
%     end

    % Inhibitory synaptic tau constants
    inh_tau1 = taui_dec;
    inh_tau2 = taui_rise*taui_dec / (taui_rise + taui_dec);
    inh_B = 1/((inh_tau2/inh_tau1)^(taui_rise/inh_tau1) - ...
        (inh_tau2/inh_tau1)^(taui_rise/inh_tau2));


    spikeRecord_e = sparse(Ne, size(times, 2));
    spikeRecord_i = sparse(Ni, size(times, 2));

%    display 'Simulation running...'
    t = 0;

    fired_e = zeros(Ne, 1);
    fired_i = zeros(Ni, 1);
    refrac_e_cnt = zeros(Ne, 1);
    refrac_i_cnt = zeros(Ni, 1);
    f_i = 1;
    t_spike = [];

    t_i = 1;
    for t = times
        fired_e = Ve > Vt_e;
        fired_i = Vi > Vt_i;

        Ve(fired_e) = Vr_e;
        Vi(fired_i) = Vr_i;
        
        %refrac_e_cnt(fired_e) = refrac_e(fired_e);
        %refrac_i_cnt(fired_i) = refrac_i(fired_i);
        g_adapt_e(fired_e) = g_adapt_e(fired_e) + o.refrac_e_g_inc;
        g_adapt_i(fired_i) = g_adapt_i(fired_i) + o.refrac_i_g_inc;

        spikeRecord_e(:, t_i) = double(fired_e);
        spikeRecord_i(:, t_i) = double(fired_i);
        
        gi_sum = Mi*fired_i*wi*inh_B;
        gi1 = gi1 + gi_sum;
        gi2 = gi2 + gi_sum;
        %gi = gi + ;  % No weight distribution for now here
        ge = ge + Me*fired_e;
        Isyn_e = (gi1 - gi2).*(V_rev_i - Ve);
        Isyn_i = ge.*(V_rev_e - Vi);

        
        Vmon_e(:, t_i) = Ve(Emon_i);
        Vmon_i(:, t_i) = Vi(Imon_i);
        Vmon_e(fired_e(Emon_i), t_i) = o.spikeVm;
        Vmon_i(fired_i(Imon_i), t_i) = o.spikeVm;
        Vmon_Isyn_e(:, t_i) = (gi1(Emon_i) - gi2(Emon_i))*(V_rev_i - Vclamp);
        Vmon_Isyn_i(:, t_i) = ge(Imon_i)*(V_rev_e - Vclamp);


        % Check if neurons fired and add to syn. conductances
        dVe = dt * ((El_e - Ve).*(1 + Rm_e*g_adapt_e) + Rm_e*Isyn_e + Rm_e*Ie) / taum_e;
        dVi = dt * ((El_i - Vi).*(1 + Rm_i*g_adapt_i) + Rm_i*Isyn_i + Rm_i*Ii) / taum_i;
        
        % Do not update states of cells which are in refractory period
        %dVe(refrac_e_cnt > 0) = 0;
        %dVi(refrac_i_cnt > 0) = 0;

        dge = dt * -1/taue * ge;
        dgi1 = dt * -1/inh_tau1 * gi1;
        dgi2 = dt * -1/inh_tau2 * gi2;

        Ve = Ve + dVe + noise_sigma*randn(Ne, 1);
        Vi = Vi + dVi + noise_sigma*randn(Ni, 1);
        ge = ge + dge;
        gi1 = gi1 + dgi1;
        gi2 = gi2 + dgi2;
        
        g_adapt_e = g_adapt_e + dt * -1./tau_adapt_e .* g_adapt_e;
        g_adapt_i = g_adapt_i + dt * -1./tau_adapt_i .* g_adapt_i;
        
        Ie = Ie + o.dIe*dt;
        Ii = Ii + o.dIi*dt;
        
        % Decrease refractory counters
        %refrac_e_cnt = refrac_e_cnt - 1;
        %refrac_i_cnt = refrac_i_cnt - 1;

        t_i = t_i + 1;
    end


    Vmon.e = Vmon_e;
    Vmon.i = Vmon_i;
    Vmon.Isyn_e = Vmon_Isyn_e;
    Vmon.Isyn_i = Vmon_Isyn_i;
    Vmon.t = Vmon_t;

end
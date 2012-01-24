function opt = getStdCellProperties()
        % Excitatory cells
        opt.taum_e = 9.3e-3;
        opt.taue = 1e-3;
        opt.El_e = -68.5e-3;
        opt.Vt_e = -50.0e-3;
        opt.Vr_e = opt.El_e;
        opt.Rm_e = 44e6; % MOhm
        opt.refrac_e_mean = 40e-3;
        opt.refrac_e_std = 5e-3;
        opt.refrac_e_g_inc = 1/opt.Rm_e/2;


        % Inhibitory cell
        opt.taum_i = 10e-3;
        opt.taui_dec = 5e-3;
        opt.taui_rise = 1e-3;
        opt.El_i = -60e-3;
        opt.Vt_i = -50e-3;
        opt.Vr_i = opt.El_i;
        opt.Rm_i = 44e6; % MOhm
        opt.refrac_i_mean = 7.5e-3; %msec
        opt.refrac_i_std  = 0.5e-3;
        opt.refrac_i_g_inc = 1/opt.Rm_i;

        opt.spikeVm = 40e-3;
        
        % Reversal potentials
        opt.V_rev_e = 0e-3;
        opt.V_rev_i = -75e-3;
        
        opt.Vclamp = -50e-3;
end
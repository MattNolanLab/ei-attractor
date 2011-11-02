% Comparison of F-I curves of integrate and fire neurons with/without
% refractory period

close all;
clear all;

refrac = (0:2:14) * 1e-3;

pA = 1e12;
fontSize = 15;


% Excitatory cells
E_cell.taum = 9.3e-3;
E_cell.El = -68.5e-3;
E_cell.Vt = -50.0e-3;
E_cell.Vr = E_cell.El;
E_cell.Rm = 44e6; % MOhm


% Inhibitory cell
I_cell.taum = 2e-3;
I_cell.El = -60e-3;
I_cell.Vt = -50e-3;
I_cell.Vr = I_cell.El;
I_cell.Rm = 44e6; % MOhm


I0 = 0;
Imax = 900e-12; %pA


figure('Position', [800 0 1000 500]);

% Stellate cell
subplot(1,2,1, 'FontSize', fontSize);
hold all;
for r_it = 1:numel(refrac)
    [T I] = spikingPeriodIF(E_cell, refrac(r_it), I0, Imax);
    F = 1./T;

    F_i = find(imag(T) ~= 0);
    F(F_i) = 0;

    plot(I*pA, F, 'LineWidth', 1.5);
    leg{r_it} = sprintf('\\Delta_r = %d ms', refrac(r_it)*1000);
end

xlabel('Current (pA)');
ylabel('Firing rate (Hz)');
legend(leg, 'Location', 'NorthWest', 'FontSize', fontSize - 2);
title('A');
ylim([0 200]);
box on;


% Interneuron
subplot(1,2,2, 'FontSize', fontSize);
hold all;
for r_it = 1:numel(refrac)
    [T I] = spikingPeriodIF(I_cell, refrac(r_it), I0, Imax);
    F = 1./T;

    F_i = find(imag(T) ~= 0);
    F(F_i) = 0;

    plot(I*pA, F, 'LineWidth', 1.5);
    leg{r_it} = sprintf('\\Delta_r = %d ms', refrac(r_it)*1000);
end

xlabel('Current (pA)');
ylabel('Firing rate (Hz)');
legend(leg, 'Location', 'NorthWest', 'FontSize', fontSize-2);
title('B');
ylim([0 1000]);
box on;


set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
print('-depsc2', 'F_I_curves_stellate_interneuron.eps');

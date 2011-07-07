% Plot the linear response of a passive chunk of integrate and fire
% membrane
close all;

tau = 0.01; % s
R = 1; % 10 MOhm

f = 0:0.01:150; % Hz

% Complex frequency response
H_f = R./(1 + i*2*pi*tau.*f);


fontSize = 14;
figure('Position', [600 500 1000 600]);

subplot(2, 1, 1, 'FontSize', fontSize);
plot(f, abs(H_f), 'k');
ylabel('|H(f)|');
title('Frequency response of a passive membrane');

subplot(2, 1, 2, 'FontSize', fontSize);
plot(f, angle(H_f)/2/pi*360, 'k');
xlabel('Frequency (Hz)');
ylabel('\theta (f)');


%
% Filtering properties given HCN+/+ and HCN-/-
%

% [HCN+/+ HCN-/-]

R_pp = 55.2 * 1e6;
R_mm = 87.3 * 1e6;
tau_pp = 10.3 * 1e-3;
tau_mm = 19.8 * 1e-3;
f = 0:0.01:150;

H_pp = R_pp./(1 + i*2*pi*tau_pp.*f);
H_mm = R_mm./(1 + i*2*pi*tau_mm.*f);

% Complex impedance
A = [H_pp; H_mm] * 1e-6;

fontSize = 14;
figure('Position', [600 500 1000 600]);


subplot(2, 1, 1, 'FontSize', fontSize, 'XMinorTick','on');
box on;
hold all;
plot(f, abs(A));
ylabel('|H(f)| (M\Omega)');
title('Frequency response of HCN1+/+ vs HCN-/- stellate cell');
legend('HCN1^{+/+}', 'HCN1^{-/-}');

subplot(2, 1, 2, 'XMinorTick','on', 'FontSize', fontSize);
box on;
hold all;
plot(f, angle(A)/2/pi*360);
xlabel('Frequency (Hz)');
ylabel('\theta (f)');
legend('HCN1^{+/+}', 'HCN1^{-/-}');


%
% Filtering properties of a resonant pasive membrane
%
tau  = 10.3 * 1e-3; % s
tau1_list = [20:20:120] * 1e-3;
R = 55.2*1e6; % 10 MOhm
gamma = 1;


fontSize = 14;
figure('Position', [600 500 1000 600]);

leg = {};
it = 1;

for tau1 = tau1_list
    f = 0:0.01:50; % Hz

    % Complex frequency response
    H_f = R./(1 + i*2*pi*tau.*f + gamma./(1 - i*2*pi.*f*tau1)) / 1e6;



    subplot(2, 1, 1, 'FontSize', fontSize);
    hold all;
    plot(f, abs(H_f));
    ylabel('|H(f)| (M\Omega)');
    title(['LIF ECII stellate cell fit with resonant variable (tau = ' num2str(tau*1e3) ' ms', ', R = ' num2str(R/1e6) ' M\Omega)']);

    subplot(2, 1, 2, 'FontSize', fontSize);
    hold all;
    plot(f, angle(H_f)/2/pi*360);
    xlabel('Frequency (Hz)');
    ylabel('\theta (f)');
    leg{it} = ['\tau_1 = ' num2str(tau1*1e3) ' ms'];
    it = it + 1;
end

subplot(2, 1, 1);
legend(leg); %, 'Location', 'EastOutside');

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'stellate_cell_resonant_filter.eps');

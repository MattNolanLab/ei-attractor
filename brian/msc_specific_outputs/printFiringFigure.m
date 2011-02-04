% Print a figure illustrating firing of a single neuron

load results/job23600_2010-08-10T16-50-33_output.mat

fontSize = 18;

SNList_nID = 10;

figure(1);
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*1000);

x_limits = [9.2 9.6];

figure(2);
subplot(1, 1, 1, 'FontSize', fontSize);
plot((SNMonitor_times - x_limits(1))*1000, SNMonitor_values(SNList_nID, :)*1000, 'k');
xlim((x_limits - x_limits(1))*1000);

xlabel('Time (ms)');
ylabel('Membrane potential (mV)');

print -depsc2 '../../thesis/src/fig/IF_firing.eps'
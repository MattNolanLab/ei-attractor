close all;
clearvars -except c00*


fontSize = 14;
figure('Position', [1200 500 1000 600]);

dt = c001_Time(2) - c001_Time(1);

% plot stellate cell
subplot(2,1,1, 'FontSize', fontSize);
plot(c001_Time, c003_Membrane_Voltage_2*1000, 'r');
box off;
ylabel('V_m (mV)');


subplot(2,1,2, 'FontSize', fontSize);
plot(c001_Time, c002_Membrane_Voltage_1*1000);
box off;
xlabel('Time (s)');
ylabel('V_m (mV)');
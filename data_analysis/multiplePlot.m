% multiple plot

for jid = 12200:12251
    d = dir(['results/fig/batch/net_stability/connMult4_lambda_s80/job' num2str(jid) '*.mat']);
    plotStatistics(['results/fig/batch/net_stability/connMult4_lambda_s80/' d(end).name], [2048]);
end
% Delete very small interactions between neurons (weights) and make the
% matrix sparse

nS = 1e-9;
threshold = 0.001*nS;

conn_sp = connections;
conn_sp(find(conn_sp < threshold)) = 0;
conn_sp = sparse(conn_sp);
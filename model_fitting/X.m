function res = X(x, y, j, k)

N = numel(x);

if (N ~= numel(y))
    error('x and y inputs must be the same length')
end

if (j < 0 || k < 0)
    error('j and k must be >= 0');
end


x = [zeros(1, j) x(1:end-j)];
y = [zeros(1, k) y(1:end-k)];

res = 1/N * sum(x .* y) - 1/N^2 * sum(x) * sum(y);

end
function res = S(x, y, j)

N = numel(x);

if (N ~= numel(y))
    error('x and y inputs must be the same length')
end

if (j < 0)
    error('j must be >= 0');
end


y = [zeros(1, j) y(1:end-j)];

res = 1/N * sum(x .* y) - 1/N^2 * sum(x) * sum(y);

end
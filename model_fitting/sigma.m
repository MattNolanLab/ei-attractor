function res = sigma(x, y)

N = numel(x);

if (N ~= numel(y))
    error('x and y inputs must be the same length')
end

res = 1/N*sum(x .* y) - 1/N^2*sum(x)*sum(y);

end
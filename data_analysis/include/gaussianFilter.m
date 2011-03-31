function [result] = gaussianFilter(X, h)
    result = exp(-X.^2/h^2);
end
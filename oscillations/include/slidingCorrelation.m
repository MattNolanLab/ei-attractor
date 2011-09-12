function [C] = slidingCorrelation(s1, s2, win_len)
    % Assuming for simplicity both signals have the same length and are row
    % vectors
    
    len = min(size(s1, 2), size(s2, 2));
    C = zeros(1, len-win_len)*nan;
    
    for it = 1:len-win_len
        tmpc = corrcoef(s1(it:it+win_len), s2(it:it+win_len));
        C(it) = tmpc(1, 2);
    end
end
function [r, c] = getBlobNeighborhood(blobPos_r, blobPos_c, rad, sheet_size)
% compute an array of neurons' 2D positions - a neighbourhood with radius
% rad

    [r, c] = meshgrid(-rad:rad);
    r = reshape(r, 1, (2*rad+1)^2) + blobPos_r;
    c = reshape(c, 1, (2*rad+1)^2) + blobPos_c;
    
    % Make the neighborhood circular
    id_del = find(sqrt((r - blobPos_r).^2 + (c - blobPos_c).^2) > rad);
    r(id_del) = [];
    c(id_del) = [];
    
    % Wrap around in x and y directions
    r_wrap = find(r < 0); r(r_wrap) = sheet_size + r(r_wrap);
    c_wrap = find(c < 0); c(c_wrap) = sheet_size + c(c_wrap);
    
    r_wrap = find(r > sheet_size - 1);
    r(r_wrap) = r(r_wrap) - sheet_size;
    
    c_wrap = find(c > sheet_size - 1);
    c(c_wrap) = c(c_wrap) - sheet_size;
end
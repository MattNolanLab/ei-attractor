function [r, c] = getBlobNeighborhood(blobPos_r, blobPos_c, rad)
% compute an array of neurons' 2D positions - a neighbourhood with radius
% rad

    [r, c] = meshgrid(-rad:rad);
    r = reshape(r, 1, (2*rad+1)^2) + blobPos_r;
    c = reshape(c, 1, (2*rad+1)^2) + blobPos_c;

end
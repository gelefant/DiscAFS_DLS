function V = chebVandInt(deg,centers,radii)
% -------------------------------------------------------------------------
% It computes the Vandermonde matrix for the interpolant constructs by 
% integrals on discs, via Chebyshev polynomials   
% 
% INPUT:
% deg     - degree of the interpolant
% centers - a matrix Nx2 of the coordinates of the centers in the unitarian
%           disc
% radii   - a column vector di dimension N of the radii of the discs
% OUTPUT
% V       - Vandermonde matrix
% -------------------------------------------------------------------------
% Dates
% --------------------------------------------------------------------------
% First version: November 15, 2023;
% Checked: December 07, 2023.
% --------------------------------------------------------------------------
% Authors
% --------------------------------------------------------------------------
% L. Bruni Bruno and G. Elefante
% --------------------------------------------------------------------------
% Paper
% --------------------------------------------------------------------------
% "Interpolation by integrals on discs"
% L. Bruni Bruno and G. Elefante
% --------------------------------------------------------------------------
dimP = nchoosek(deg+2,2);

[XYW,dbox]=cub_unit_disk_sets(deg);

PolDeg = polydeg(deg); 

X = XYW(:,1); Y = XYW(:,2); W = XYW(:,3);

for i = 1:length(radii)
    VxCheb = chebpolys(deg,X*radii(i)+centers(i,1));
    VyCheb = chebpolys(deg,Y*radii(i)+centers(i,2));
    for j = 1:dimP
        V(i,j) = VxCheb(:,PolDeg(j,1)+1)'.*VyCheb(:,PolDeg(j,2)+1)'*W*radii(i)^2;
    end
end

end
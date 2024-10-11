function V = chebVand(deg,centers)
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
%--------------------------------------------------------------------------
% First version: November 15, 2023;
% Checked: December 07, 2023.
%--------------------------------------------------------------------------
% Authors
%--------------------------------------------------------------------------
% L. Bruni Bruno and G. Elefante
%--------------------------------------------------------------------------
% Paper
%--------------------------------------------------------------------------
% "Interpolation by integrals on discs"
% L. Bruni Bruno and G. Elefante
%--------------------------------------------------------------------------
dimP = nchoosek(deg+2,2);

PolDeg = polydeg(deg); 

for i = 1:size(centers,1)
    VxCheb = chebpolys(deg,centers(i,1));
    VyCheb = chebpolys(deg,centers(i,2));
    for j = 1:dimP
        V(i,j) = VxCheb(:,PolDeg(j,1)+1)'.*VyCheb(:,PolDeg(j,2)+1);
    end
end

end
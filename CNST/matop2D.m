%=========================================================================
% Copyright (C) University of Santiago de Compostela
% Author: Saray Busto Ulloa
%=========================================================================
% Licensed under the Creative Commons Attribution-NonCommercial 4.0
% International licence (CC BY-NC 4.0):
% https://creativecommons.org/licenses/by-nc/4.0/
% You may use, share and adapt this file for non-commercial purposes only,
% provided that you credit the authors and the source, and indicate any
% changes you make.
%=========================================================================
% Funzione per il calcolo del prodotto matrice vettore per
% risolvere il sistema di pressioni delle equazioni di Navier-Stokes inc. 2D
function  Ap = matop2D(p)
global dx dy IMAX JMAX
Ap = zeros(IMAX,JMAX);
for i=1:IMAX
    for j=1:JMAX
        % x-direction
        if (i==1)
            Ap(i,j) = ( (p(i+1,j)-p(i,j))/dx - p(i,j)/(dx/2) )/dx;
        elseif (i==IMAX)
            Ap(i,j) = ( -p(i,j)/(dx/2) - (p(i,j)-p(i-1,j))/dx )/dx;
        else
            Ap(i,j) = ( (p(i+1,j)-p(i,j))/dx - (p(i,j)-p(i-1,j))/dx )/dx;
        end
        % y-direction
         if (j==1)
             Ap(i,j) = Ap(i,j) + ( (p(i,j+1)-p(i,j))/dy - p(i,j)/(dy/2) )/dy;
        elseif (j==JMAX)
            Ap(i,j) = Ap(i,j) + ( -p(i,j)/(dy/2) - (p(i,j)-p(i,j-1))/dy )/dy;
        else
            Ap(i,j) = Ap(i,j) + ( (p(i,j+1)-p(i,j))/dy - (p(i,j)-p(i,j-1))/dy )/dy;
        end
    end
end

end
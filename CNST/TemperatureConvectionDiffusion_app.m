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
%% Funzione per il calcolo dei termini convettivi
%    e diffusivi per la equazione della temepratura
%
% Copy of TemperatureConvectionDiffusion.m used by CNS.mlapp, with two extra
% boundary conditions (everything else is identical):
%   - Tprofile = 'Adiabática': adiabatic bottom wall (no heat flux);
%   - radMask/Trad: fixed temperature Trad on the cells of the left wall where
%     radMask is true (idealised radiator on the lower third of the wall);
%   - winMask/Twin: fixed temperature Twin on the cells of the left wall where
%     winMask is true (idealised window on the middle third of the wall).
function Tnew= TemperatureConvectionDiffusion_app(T,u,v,kappa)
global dx dy IMAX JMAX dt Tb xb Tprofile radMask Trad winMask Twin
Tnew = T;
for i=1:IMAX
    for j=1:JMAX
        % x flux
        if (i==1)
            fp = 0.5*u(i+1,j)*(T(i,j)+T(i+1,j))  - 0.5*abs(u(i+1,j))*(T(i+1,j)-T(i,j)) - kappa *(T(i+1,j)-T(i,j))/dx;
            if radMask(j)
                fm = -kappa*(T(i,j)-Trad)/(dx/2); % Fixed temperature on the radiator
            elseif winMask(j)
                fm = -kappa*(T(i,j)-Twin)/(dx/2); % Fixed temperature on the window
            else
                fm = 0; % adiabatic wall (non c'è flusso di calore in questo boundary)
            end
        elseif (i==IMAX)
            fp = 0; % adiabatic wall
            fm = 0.5*u(i,j)*(T(i,j)+T(i-1,j)) - 0.5*abs(u(i,j))*(T(i,j)-T(i-1,j)) - kappa*(T(i,j)-T(i-1,j))/dx;
        else
            fp = 0.5*u(i+1,j)*(T(i,j)+T(i+1,j))  - 0.5*abs(u(i+1,j))*(T(i+1,j)-T(i,j)) - kappa *(T(i+1,j)-T(i,j))/dx;
            fm = 0.5*u(i,j)*(T(i,j)+T(i-1,j)) - 0.5*abs(u(i,j))*(T(i,j)-T(i-1,j)) - kappa*(T(i,j)-T(i-1,j))/dx;
        end
        % y flux
        if (j==1)
            gp = 0.5*v(i,j+1)*(T(i,j+1)+T(i,j)) - 0.5*abs(v(i,j+1))*(T(i,j+1)-T(i,j)) - kappa*( T(i,j+1)-T(i,j))/dy;
            if strcmp(Tprofile,'Adiabática')
                gm = 0; % adiabatic bottom wall
            else
                if strcmp(Tprofile,'Lineal')
                    TBC = Tb ;%*exp(-0.5*(xb(i)-0.5)^2/0.1^2); % Fixed temperature at the bottom
                else
                    TBC = Tb*exp(-0.5*(xb(i)-0.5)^2/0.1^2); % Fixed temperature at the bottom
                end
                gm = -kappa*(T(i,j)-TBC)/(dy/2);
            end
        elseif (j==JMAX)
            gp = 0;
            gm = 0.5*v(i,j)*(T(i,j)+T(i,j-1)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,j-1)) - kappa*(T(i,j)-T(i,j-1))/dy;
        else
            gp = 0.5*v(i,j+1)*(T(i,j+1)+T(i,j)) - 0.5*abs(v(i,j+1))*(T(i,j+1)-T(i,j)) - kappa*( T(i,j+1)-T(i,j))/dy;
            gm = 0.5*v(i,j)*(T(i,j)+T(i,j-1)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,j-1)) - kappa*(T(i,j)-T(i,j-1))/dy;
        end
        Tnew(i,j) = Tnew(i,j) - dt/dx*(fp-fm) -dt/dy*(gp-gm);
    end
end
end

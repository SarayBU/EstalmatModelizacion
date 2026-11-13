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
% Function for the computation of convective and viscous terms 
% for the incompressible Navier-Stokes equations
%
% Copy of MomentumConvectionDiffusion.m used by CNS.mlapp. The only change is
% that the Boussinesq buoyancy term uses the reference temperature Tref
% (the initial temperature): beta*g*(T - Tref). With Tref = 0 it is identical
% to the original.
function [ustar,vstar] = MomentumConvectionDiffusion_app(u,v,T)
global dx dy IMAX JMAX dt  g  nu  beta uLid Tref
ustar= u;
vstar = v;
%% Compute ustar
for i=1:IMAX+1
    for j=1:JMAX
        % x direction (u*u_x+nu*u_xx)
        if (i==1)
            ustar(i,j) = 0; % no flow through the wall
        elseif (i==IMAX+1)
            ustar(i,j) = 0; % no flow through the wall
        else
            % Convective terms (u*u_x)
            uav = 0.5*(u(i,j)+u(i-1,j)); % Averaged velocity
            um = 0.5*(uav-abs(uav)); % Left-going waves
            uav = 0.5*(u(i+1,j)+u(i,j)); % Averaged velocity
            up = 0.5*(uav+abs(uav)); % Right-going waves
            ustar(i,j) = ustar(i,j) -dt/dx*(um*(u(i+1,j)-u(i,j)) + up*(u(i,j)-u(i-1,j)) );
            % Viscous terms (nu*u_xx)
            ustar(i,j) = ustar(i,j) + nu*dt/dx*( (u(i+1,j)-u(i,j))/dx - (u(i,j)-u(i-1,j))/dx );
        end        
        % y direction   (v*u_y+nu*u_yy) 
        if (j==1)
            % Convective terms (v*u_y)
            vav = 0.5*(v(min(i,IMAX),j+1)+v(max(1,i-1),j+1)); % Average velocity
            vm = 0.5*(vav-abs(vav)); % Downward moving waves
            ustar(i,j) = ustar(i,j) -dt/dy*(vm*(u(i,j+1)-u(i,j)) + 0 ); % no flux from below
            % Viscous terms (nu*u_yy)
            ustar(i,j) = ustar(i,j) + nu*dt/dy*( (u(i,j+1)-u(i,j))/dy - (u(i,j)-0)/(dy/2) ); % no-slip bC
        elseif (j==JMAX)
            % Convective terms (v*u_y)
            vav = 0.5*(v(min(i,IMAX),j)+v(max(1,i-1),j)); % Average velocity
            vp = 0.5*(vav+abs(vav)); % Upward moving waves
            ustar(i,j) = ustar(i,j) -dt/dy*(0 + vp*(u(i,j)-u(i,j-1)) );
            % Viscous terms (nu*u_yy)
            ustar(i,j) = ustar(i,j) + nu*dt/dy*( (uLid-u(i,j))/(dy/2) - (u(i,j)-u(i,j-1))/dy );        
        else
            % Convective terms (v*u_y)
            vav = 0.5*(v(min(i,IMAX),j+1)+v(max(1,i-1),j+1)); % Average velocity
            vm = 0.5*(vav-abs(vav)); % Downward moving waves
            vav = 0.5*(v(min(i,IMAX),j)+v(max(1,i-1),j)); % Average velocity
            vp = 0.5*(vav+abs(vav)); % Upward moving waves
            ustar(i,j) = ustar(i,j) -dt/dy*(vm*(u(i,j+1)-u(i,j)) + vp*(u(i,j)-u(i,j-1)) );
            % Viscous terms (nu*u_yy)
            ustar(i,j) = ustar(i,j) + nu*dt/dy*( (u(i,j+1)-u(i,j))/dy - (u(i,j)-u(i,j-1))/dy );
        end
    end
end


%% Compute vstar
for i=1:IMAX
    for j=1:JMAX+1
        vstar(i,j) = vstar(i,j) + dt*beta*g*0.5*(T(i,min(JMAX,j))+T(i,max(1,j-1))) - dt*beta*g*Tref; % Boussinesq assumption (reference temperature Tref)
        % x direction (u*v_x+nu*v_xx)
        if (i==1)
            uav = 0.5*(u(i+1,min(j,JMAX))+u(i+1,max(1,j-1))); % Average velocity
            um = 0.5*(uav-abs(uav)); 
            vstar(i,j) = vstar(i,j) -dt/dx*(um*(v(i+1,j)-v(i,j)) + 0 );
            % Viscous terms (nu*v_xx)
            vstar(i,j) = vstar(i,j) + nu*dt/dx*( (v(i+1,j)-v(i,j))/dx - (v(i,j)-0)/(dx/2) );
        elseif (i==IMAX)
            uav = 0.5*(u(i,min(j,JMAX))+u(i,max(1,j-1))); % Average velocity
            up = 0.5*(uav+abs(uav)); 
            vstar(i,j) = vstar(i,j) -dt/dx*(0 + up*(v(i,j)-v(i-1,j)) );
            % Viscous terms (nu*v_xx)
            vstar(i,j) = vstar(i,j) + nu*dt/dx*( (0-v(i,j))/(dx/2) - (v(i,j)-v(i-1,j))/dx );
        else
            % Convective terms (u*v_x)
            uav = 0.5*(u(i+1,min(j,JMAX))+u(i+1,max(1,j-1))); % Average velocity
            um = 0.5*(uav-abs(uav)); % Left-going waves
            uav = 0.5*(u(i,min(j,JMAX))+u(i,max(1,j-1))); % Average velocity
            up = 0.5*(uav+abs(uav)); % Right-going waves
            vstar(i,j) = vstar(i,j) -dt/dx*(um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) );
            % Viscous terms (nu*v_xx)
            vstar(i,j) = vstar(i,j) + nu*dt/dx*( (v(i+1,j)-v(i,j))/dx - (v(i,j)-v(i-1,j))/dx );
        end            
        % y direction   (v*v_y+nu*v_yy)
        if (j==1)
            vstar(i,j) = 0; % no flow through the wall
        elseif (j==JMAX+1)            
            vstar(i,j) = 0; % no flow through the wall        
        else
            % Convective terms (v*u_y)
            vav = 0.5*(v(i,j)+v(i,j+1)); % Average velocity
            vm = 0.5*(vav-abs(vav)); % Downward moving waves
            vav = 0.5*(v(i,j-1)+v(i,j)); % Average velocity
            vp = 0.5*(vav+abs(vav)); % Upward moving waves
            vstar(i,j) = vstar(i,j) -dt/dy*(vm*(v(i,j+1)-v(i,j)) + vp*(v(i,j)-v(i,j-1)) );
            % Viscous terms (nu*u_yy)
            vstar(i,j) = vstar(i,j) + nu*dt/dy*( (v(i,j+1)-v(i,j))/dy - (v(i,j)-v(i,j-1))/dy );
        end
    end
end

end
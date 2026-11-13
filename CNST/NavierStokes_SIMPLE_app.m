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
function result = NavierStokes_SIMPLE_app(xR,yR,Nx,Ny,tend,CFL,uL,Tt,Tprof,gravity,mu,kpp,bt,radOn,Trad_in,winOn,Twin_in,T0,app)
%
% Copy of NavierStokes_SIMPLE.m used by CNS.mlapp.
% The numerical method is identical to NavierStokes_SIMPLE.m. The only changes
% are that plotting is delegated to the app (app.newFrame), which can also
% pause or stop the simulation between time steps.
% New boundary conditions (TemperatureConvectionDiffusion_app.m):
%   - Tprof = 'Adiabática': adiabatic bottom wall;
%   - radOn = true: fixed temperature Trad_in on the lower third of the left
%     wall (idealised radiator);
%   - winOn = true: fixed temperature Twin_in on the middle third of the left
%     wall (idealised window).
% New initial condition: uniform temperature T0 (0 in NavierStokes_SIMPLE.m),
% also used as reference temperature of the buoyancy term
% (MomentumConvectionDiffusion_app.m).

%% SIMPLE algorithm per cercare la soluzione delle equazioni
%    di Navier-Stokes incomprimibili in 2D (con effeti di galegiamento)
global dx dy IMAX JMAX dt g nu beta uLid xb Tb Tprofile radMask Trad winMask Twin Tref

%% Definizione del problema fisico
nu = mu;      % viscosita cinematica
kappa = kpp; % Coefficente di conduzione di calore (kappa=lambda/(rho*c))
beta = bt;           % coefficente di espansione termica
g = gravity;           % gravità
uLid = uL;          % velocità nel bordo superiore
Tb = Tt;             % temperatura al bordo inferiore
Tprofile = Tprof;
Tref = T0;          % reference temperature for the Boussinesq buoyancy

%% Definire il dominio computazionale
xL = 0; yL = 0;
% Discretizazione
IMAX = Nx; JMAX = Ny;
dx = (xR-xL)/IMAX; dy = (yR-yL)/JMAX;
% Baricentri (pressione)
xb = linspace(xL+dx/2,xR-dx/2,IMAX);
yb = linspace(yL+dy/2,yR-dy/2,JMAX);
% Vertici (velocita)
x = linspace(xL,xR,IMAX+1); % coordinate x della u
y = linspace(yL,yR,JMAX+1); % coordinate y della v
% Radiator: cells of the left wall whose centre lies in the lower third
Trad = Trad_in;
radMask = radOn & (yb <= yL + (yR-yL)/3);
% Window: cells of the left wall whose centre lies in the middle third
Twin = Twin_in;
winMask = winOn & (yb > yL + (yR-yL)/3) & (yb <= yL + 2*(yR-yL)/3);

%% Definizione delle variabili per il tempo
NMAX = 10000; % Numero massimo di iterazioni temporali
time = 0;            % Tempo iniziale
CFLt = CFL;          % Condizione CFL (per il calco del dt, per essere sicuri della stabilità)

%% Definire le condizioni iniziali
u = zeros(IMAX+1,JMAX); v = zeros(IMAX,JMAX+1); % campo di velocita
p = zeros(IMAX,JMAX); % presione
T = T0*ones(IMAX,JMAX); % temperatura (uniform initial value T0)

% Definizione della velocita nei baricentri per potter fare il plot
for i=1:IMAX
    for j=1:JMAX
        ub(i,j) = 0.5*(u(i,j)+u(i+1,j));
        vb(i,j) = 0.5*(v(i,j)+v(i,j+1));
    end
end
% Plot delle condizioni iniziali (delegated to the app)
setup = struct('xb',xb,'yb',yb,'xL',xL,'xR',xR,'yL',yL,'yR',yR,'tend',tend,'Tb',Tb, ...
    'Tprofile',Tprofile,'radOn',radOn,'Trad',Trad,'radTop',yL+dy*sum(radMask), ...
    'T0',T0,'winOn',winOn,'Twin',Twin,'winBot',yL+dy*(find(winMask,1)-1),'winTop',yL+dy*find(winMask,1,'last'));
Etot = sum(sum(T))*dx*dy;
stopped = app.newFrame(setup,time,T,ub,vb,Etot);

for n=1:NMAX % loop in tempo
    if stopped
        break
    end
    % Calcolare il passo temporale
    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));
    dt      = CFLt/( umax/dx + vmax/dy + 2*(nu+kappa)*(1/dx^2+1/dy^2) );
    % Controlliamo di arrivare al tempo finale
    if (time+dt>tend)
        dt = tend-time;
    end
    % Stop criterium
    if ( abs(time-tend) <= 1.e-10 )
        break
    end

    %% 1.- Calcolare il contributo dei termi convettivi e diffusivi per la equazione dei momenti
    [ustar,vstar] = MomentumConvectionDiffusion_app(u,v,T);

    % Discretizare la equazione della energia
    T = TemperatureConvectionDiffusion_app(T,u,v,kappa);
    Etot = sum(sum(T))*dx*dy;
    fprintf('Time = %f, Etot = %e \n',time,Etot)

    % 1.2.- Aggiungere la contribuzione della pressione
    for i=1:IMAX+1
        for j=1:JMAX
            if (i==1)
                ustar(i,j) = 0; % Impose no flow through the wall
            elseif(i==IMAX+1)
                ustar(i,j) = 0; % Impose no flow through the wall
            else
                ustar(i,j) = ustar(i,j) - dt/dx*(p(i,j)-p(i-1,j));
            end
        end
    end
    for i=1:IMAX
        for j=1:JMAX+1
            if (j==1)
                vstar(i,j) = 0; % Impose no flow through the wall
            elseif(j==JMAX+1)
                vstar(i,j) = 0; % Impose no flow through the wall
            else
                vstar(i,j) = vstar(i,j) - dt/dy*(p(i,j)-p(i,j-1));
            end
        end
    end

    %% 2.- Calcolo del rhs per il problema di Poisson per la pressione
    for i=1:IMAX
        for j=1:JMAX
            rhs(i,j) = (ustar(i+1,j)-ustar(i,j))/(dt*dx) + (vstar(i,j+1)-vstar(i,j))/(dt*dy);
        end
    end
    % Risolvere il problema di Poisson
    dp = CGop2D(rhs);

    %% 3.- Aggiornare velocita e pressione
    for i=1:IMAX+1
        for j=1:JMAX
            if (i==1)
                u(i,j) = 0;
            elseif (i==IMAX+1)
                u(i,j) = 0;
            else
                u(i,j) = ustar(i,j) - dt/dx*(dp(i,j)-dp(i-1,j));
            end
        end
    end
    for i=1:IMAX
        for j=1:JMAX+1
            if (j==1)
                v(i,j) = 0;
            elseif (j==JMAX+1)
                v(i,j) = 0;
            else
                v(i,j) = vstar(i,j) - dt/dy*(dp(i,j)-dp(i,j-1));
            end
        end
    end

    % Aggiornare la pressione
    p = p + dp;

    % Aggiornare il tempo
    time = time + dt;

    % Plot
    for i=1:IMAX
        for j=1:JMAX
            ub(i,j) = 0.5*(u(i,j)+u(i+1,j));
            vb(i,j) = 0.5*(v(i,j)+v(i,j+1));
        end
    end
    stopped = app.newFrame([],time,T,ub,vb,Etot);
end

result = struct('time',time,'completed',abs(time-tend) <= 1.e-10);
end

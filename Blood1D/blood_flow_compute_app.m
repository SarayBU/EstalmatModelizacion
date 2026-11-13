%=========================================================================
% Copyright (C) University of Santiago de Compostela
% Copyright (C) University of A Coruña
% Authors: Saray Busto Ulloa, Elena Vázquez Cendón, Andrés Prieto Aneiros
%=========================================================================
% Licensed under the Creative Commons Attribution-NonCommercial 4.0
% International licence (CC BY-NC 4.0):
% https://creativecommons.org/licenses/by-nc/4.0/
% You may use, share and adapt this file for non-commercial purposes only,
% provided that you credit the authors and the source, and indicate any
% changes you make.
%=========================================================================
function result = blood_flow_compute_app(aL_aux,aR_aux,uL,uR,NCELLS,CFL,tEnd,E,h,r0,app)
% BLOOD_FLOW_COMPUTE_APP  Copy of blood_flow_compute.m used by BloodFlow.mlapp.
%
% The numerical method (vessel, model, time loop, exact Riemann solution and
% error computation) is identical to blood_flow_compute.m. The only changes are:
%   - E, h, r0 and tEnd are inputs (defaults in the app: the former fixed values);
%   - plotting is delegated to the app (app.newFrame), which can also pause or
%     stop the simulation between time steps;
%   - the errors are returned to the app as well as printed.

%% DEFINE PARAMETERS FOR VESSELSdefine parameters for vessels
L=0.2; %m
%E=4.0e5; %pascal
%h=1.1e-3; %m
%r0=9.99e-3; %m
a0=pi*r0^2; %m^2
rho=1050; %kg/m^3
gate=0.1; %m
% Arterias
m=0.5;
K=sqrt(pi)/(1-0.5^2)*E*h/sqrt(a0);
gamma=m*K/rho/(m+1)/a0^m;

%% DEFINE THE VESSELS
v=vessel(NCELLS,L,2,K,gamma,a0,1,1,10);
vessEx = vessel(NCELLS,L,2,K,gamma,a0,1,1,10);

%% DEFINE THE MODEL
mod=model(2,CFL);
mod=mod.add(v);
mod=mod.add(vessEx);

%% INITIAL CONDITION
%tEnd = 0.014;
gate = 0.1;
aL = aL_aux*v.a0;
aR = aR_aux*v.a0;

%timing
time=0.;
dt=0.;
nMax=100000; %maximum number of time iterations

%% EXACT VESSEL
[aS,uS] = mod.solveERP(1,aL,aR,uL,uR);
exact = zeros(mod.C(2).NCELLS,3);
exact(:,1) = mod.C(2).xC';
mod.C(2).Q = mod.exactSampleERP(aL,aR,uL,uR,2,aS,uS,time,gate);
exact(:,2) = (mod.C(2).Q(1,:)/mod.C(2).a0)';
exact(:,3) = mod.C(2).Q(2,:)';

%% INITIAL CONDITION
mod.C(1).Q(1,:)=(mod.C(1).xC<=gate)*aL+(mod.C(1).xC>gate)*aR;
mod.C(1).Q(2,:)=(mod.C(1).xC<=gate)*aL*uL+(mod.C(1).xC>gate)*aR*uR;

%% FIGURE
xplot = mod.C(1).xC/mod.C(1).L;
aplot = mod.C(1).Q(1,:)'/mod.C(1).a0;
uplot = mod.C(1).Q(2,:)'./mod.C(1).Q(1,:)';

%% EIGENVALUE
cS=sqrt(mod.C(1).gamma*3/2*sqrt(aS));
cL=sqrt(mod.C(1).gamma*3/2*sqrt(aL));
cR=sqrt(mod.C(1).gamma*3/2*sqrt(aR));
T=linspace(0,tEnd,NCELLS);
sL=uL-mK(mod,gamma,aL,aS)/aL;
sR=uR+mK(mod,gamma,aR,aS)/aR;
uLplot =T.*(sL);
uRplot =T.*(sR);
uCplot =T.*(uS);
Uplot(1,:)=uLplot;
Uplot(2,:)=uRplot;
if (aS<=aR)
    [a,b]=size(Uplot);
    uRplot =T.*(uR+cR);
    uR2plot=T.*(uS+cS);
    Uplot(a+1,:)=uR2plot;
end
if (aS<=aL)
    [a,b]=size(Uplot);
    uLplot =T.*(uS-cS);
    uL2plot=T.*(uL-cL);
    Uplot(a+1,:)=uL2plot;
end
Uplot(1,:)=uLplot;
Uplot(2,:)=uRplot;

%% PLOT (delegated to the app)
setup = struct('xplot',xplot,'xexact',exact(:,1)/mod.C(2).L, ...
    'charX',(Uplot+gate)./mod.C(2).L,'charT',T,'tEnd',tEnd, ...
    'lmax',max(abs(sqrt([aL_aux,aR_aux]/pi))));
stopped = app.newFrame(setup,time,aplot,uplot,exact(:,2),exact(:,3));

for k=1:nMax
    if time>=tEnd || stopped
        break
    end

    dt=mod.timeStep(1);             %compute the time step
    dx=mod.C(1).dx;

    if time+dt>tEnd                 %ensure exact exit time matching
        dt=tEnd-time;
    end
    %evolve the solution
    [mod.C(1).Qbcs,mod.C(1).BCL,mod.C(1).BCR] = mod.boundaryConditions(1,mod.C(1).Q,mod.C(1).NCELLS,mod.C(1).BCL,mod.C(1).BCR,mod.C(1).Qbcs,1/2*(aL+aR),1/2*(aL+aR));
    mod.C(1).Q = mod.evolve(1,dt,dx,mod.C(1).gamma);
    %update the time
    time=time+dt;

    mod.C(2).Q = mod.exactSampleERP(aL,aR,uL,uR,2,aS,uS,time,gate);
    exact(:,2) = (mod.C(2).Q(1,:)/mod.C(2).a0)';
    exact(:,3) = mod.C(2).Q(2,:)';
    % Refresh plot
    aplot = mod.C(1).Q(1,:)'/mod.C(1).a0;
    uplot = mod.C(1).Q(2,:)'./mod.C(1).Q(1,:)';

    stopped = app.newFrame([],time,aplot,uplot,exact(:,2),exact(:,3));
end
error=zeros(2,NCELLS);
error(1,:)=abs(mod.C(2).Q(1,:)-mod.C(1).Q(1,:));
error(2,:)=abs(mod.C(2).Q(1,:).*mod.C(2).Q(2,:)-mod.C(1).Q(2,:));
l1error=zeros(2,1);
l1error(1)=mod.C(1).dx*sum(error(1,:));
l1error(2)=mod.C(1).dx*sum(error(2,:));
disp(['l1error q_1=A =',num2str(l1error(1))])
disp(['l1error q_2=Au=',num2str(l1error(2))])
%
linferror=zeros(2,1);
linferror(1)=max(error(1,:));
linferror(2)=max(error(2,:));
disp(['linferror(1)=',num2str(linferror(1))])
disp(['linferror(2)=',num2str(linferror(2))])

result = struct('l1error',l1error,'linferror',linferror,'time',time, ...
    'completed',time>=tEnd);
end

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
% Funzione per il metodo del matrix free conjugate gradient 
% per risolvere il sistema lineare A*x=b (2D, A symmetric and positive definite)
function x=CGop2D(b)
N = numel(b);               % Get the size of the system
x = b;                             % Initial guess
r = b-matop2D(x);        % Initial residual== steepest descendent direction 
p = r;                             % Initial search direction == steepest descendent direction
alpha = sum(sum(r.*r));  % square of residual
alpha0 = alpha;
tol = 1e-6;                     % Tolerance
if (alpha<tol)
    return                        % The initial guess is the exact solution
end
for k=1:N
    if (sqrt(alpha)<tol*sqrt(alpha0))
        % We have reached the expected relative tolerance
        return
    end
    v= matop2D(p);                        % lambda minimizing the 1D problem
    lambda = alpha/sum(sum(p.*v));
    x = x+lambda*p;                      % update solution
    r = r-lambda*v;                        % update residual
    alphanew = sum(sum(r.*r));
    p = r + alphanew/alpha*p;     % new search direction
    alpha = alphanew;    
end
res = sqrt(alpha);
fprintf('CG did not converge. res = %f',res)
end
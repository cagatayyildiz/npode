function [f,g] = np_ode_target(gp,par)
% NP_ODE_TARGET Computes the value and the gradients of the negavite log
% posterior, which is the target function that is minimized.
%
% INPUT
%       gp 
%       pars - optional parameters
%
% OUTPUT
%       f - negative log-posterior
%       g - gradients of f

if ~exist('par','var')
    par = par2vec(gp);
end
[f,g,~] = np_ode_fg(gp,par);
f = -f;
g = -g;
end
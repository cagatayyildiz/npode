function X = np_ode_mean_path(gp,t,x0)
% NP_ODE_MEAN_PATH Computes the mean path
%
% INPUT
%       gp 
%       t - time indices (cell array)
%       x0 - initial value (row vector)
%
% OUTPUT
%       Xs - mean path

X = np_ode_sens(t,x0,gp.F,gp.X,gp.ell,gp.sf,gp.tol);
end


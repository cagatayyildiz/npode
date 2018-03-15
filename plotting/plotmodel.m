function [] = plotmodel(gp,t)
% PLOTMODEL Plots the vector field and infered trajectories
%
% If the time points are given as input, then the ODE system is evaluated
% at these points. If there is a single trajectory to be predicted, t
% should be a column vector. In case of multiple trajectories, t is a cell
% array of column vectors.
%
% INPUT
%       gp  - npODE instance storing underlying parameters
%       t   - (optional) time points to be evaluated

if exist('t','var') 
    if ~iscell(t) && numel(t)==length(t) 
        t = {t(:)};
    end
    gp.t = t;
    gp.Nt = length(gp.t);
end

if gp.D == 2
    subplot(2,2,[1 3]);
    plotvf(gp,1);
    subplot(2,2,2);
    plotode(gp,1);
    if gp.Nt > 1
        subplot(2,2,4);
        plotode(gp,2);
    end
else
    plotode(gp,1);
end
    
end



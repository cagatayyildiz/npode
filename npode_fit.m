function gp = npode_fit(t,Y,varargin)
% NPODE_FIT Implements npODE algorithm
% The function requires two arguments Y and t (observations and associated
% time points). The rest are optional arguments, see below for detail.
%
% If there is a single input trajectory, t should be a column vector. 
% In case of multiple trajectories t is a cell array of column vectors.
% Similarly, Y is either a matrix whose rows store the observations or it
% is a cell array of matrices.
%
% INPUT 
%       t       - observation times (column or cell array)
%       Y       - observations (matrix of cell array)
%       optpars - parameters to be optimized. could be any concatenation of
%                 'log_sn', 'log_sf', 'x0', 'Fw'
%       W       - width of the inducing point grid
%       sf      - signal variance of the kernel
%       ell     - length scale of the kernel
%       model   - true ode system ('vdp', 'fhn' or 'lv')
%       x0_true - true initial value(s)
%
% OUTPUT
%       D       - problem dimensionality - int
%       x0      - initial values - Nt x 2 matrix
%       t       - time points in each trajectory - cell array of vectors
%       X       - latent states - cell array of matrices, rows store states
%       Y       - noisy obs. - cell array of matrices, rows store obs
%       odefun  - ode function that has generated the data

addpath(genpath('.'))
if ~exist('varargin','var')
    varargin = [];
end
if ~iscell(t) && numel(t)==length(t) 
    t = {t(:)};
    Y = {Y};
end

%% default options
D = size(Y{1},2);
opts.optpars    = 'log_sn-x0-Fw';
opts.W          = 6;
opts.sf         = 1;
opts.ell        = 2 * ones(1,D);
for v=1:length(varargin)/2
    opts.(varargin{v*2-1}) = varargin{v*2};
end

%% init np_ode model
npode_opts.kernel     = 'gauss'; % we only implemented sq. exp. kernel
npode_opts.Xlocs      = 'grid'; % inducing point locations
npode_opts.optpars    = opts.optpars;
npode_opts.W          = opts.W;
npode_opts.log_ell    = log(opts.ell);
npode_opts.log_sf     = log(opts.sf);

gp = np_de_model(t,Y,npode_opts);
fprintf('Initial log-posterior: %f\n',np_ode_fg(gp));

%% grid search for pre-initialization
pars = gridsearch(gp,@np_ode_fg);
gp = vec2par(gp,pars);

%% optimization
[pars,gp] = optim_npde(gp,@np_ode_target,100);
fprintf('Log-posterior after optimization: %f\n', np_ode_fg(gp));

end


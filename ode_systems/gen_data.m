function [t,Y,D,x0,X,odefun] = gen_data(varargin)
% GEN_DATA Generates noisy data using an oscillator
% The function does not require any argument but see the documentation for
% possible inputs
%
% INPUT 
%       model   - data generating oscillator - 'vdp', 'fhn' or 'lv'
%       Nt      - number of trajectories to be generated 
%       Ny      - number of data points in each trajectory
%       sn      - std of the white noise
%       t       - the duration of each trajectory. could be an integer or a
%                 vector of length Nt
%       x0      - initial values of each trajectory - of size Nt x 2
%
% OUTPUT
%       t       - time points in each trajectory - cell array of vectors
%       Y       - noisy obs. - cell array of matrices, rows store obs
%       D       - problem dimensionality - int
%       x0      - initial values - Nt x 2 matrix
%       X       - latent states - cell array of matrices, rows store states
%       odefun  - ode function that has generated the data

% default settings
opts.model  = 'vdp';
opts.Nt     = 3;
opts.Ny     = 25;
opts.sn     = 0.1;
opts.x0     = [];
D           = 2;

for v=1:length(varargin)/2
    opts.(varargin{v*2-1}) = varargin{v*2};
end

%% set properties that are specific to each time series
% t     - length of the trajectory
% Ny    - number of data points in each trajectory
% x0    - initial values
if ~isfield(opts,'t')
    if strcmp(opts.model,'vdp')
        opts.t = 7;
        odefun = @(t,x) ode_vdp(t,x);
    elseif strcmp(opts.model,'fhn')
        opts.t = 10;
        odefun = @(t,x) ode_fhn(t,x);
    elseif strcmp(opts.model,'lv')
        odefun = @(t,x) ode_lv(t,x);
        opts.t = 5;
    end
end
if ~isempty(opts.x0)
    opts.Nt = size(opts.x0,1);
end
if isscalar(opts.Ny)
    opts.Ny = opts.Ny * ones(opts.Nt,1);
end
if isscalar(opts.t)
    opts.t = opts.t * ones(opts.Nt,1);
end

% initial points
if ~isfield(opts,'x0') || isempty(opts.x0)
    tmp = lhsdesign(opts.Nt,D);
    if strcmp(opts.model,'vdp')
        x0 = tmp .* [4,4] - [2,2];
    elseif strcmp(opts.model,'fhn')
        x0 = tmp .* [6,4] - [3,2];
    elseif strcmp(opts.model,'lv')
        x0 = tmp .* [4,4] + [6,1];
    end
else
    x0 = opts.x0;
end

%% generate data
X = cell(opts.Nt,1);
Y = cell(opts.Nt,1);
t = cell(opts.Nt,1);

for i = 1:opts.Nt
    t{i} = linspace(1e-3,opts.t(i),opts.Ny(i))';
    [~,X_i] = ode45(odefun, t{i}, x0(i,:));
    X{i} = X_i;
    Y{i} = X_i + opts.sn*randn(size(X_i));
end

end
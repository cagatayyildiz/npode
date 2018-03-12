%% 
clear
close all
addpath(genpath('.'))

%% data generation
rng(5436546)
data_opts.model = 'vdp'; % oscillators. could be 'vdp', 'fhn', 'lv'
data_opts.Nt = 3; % number of input time series
data_opts.Ny = 25; % number of data points in input time series
data_opts.sn_true = 0.1; % true noise added to the trajectories

[D,x0,t,~,Y,ode_fun] = gen_data(data_opts);

%% init np_ode model
opts.x0_true    = x0;
opts.kernel     = 'gauss';
opts.optpars    = 'log_sn-x0-Fw'; % optimized parameters
opts.ode_fun    = ode_fun;
opts.Xlocs      = 'grid';  % 'minbb' 'traj' 'grid'
opts.W          = 5; % width of the inducing point grid
opts.ell        = [1.5 1.8];
opts.sf         = 0.9;

gp = np_de_model(t,Y,opts);
fprintf('Initial posterior: %f\n',np_ode_fg(gp));
par = par2vec(gp);

figure('units','normalized','outerposition',[0 0 1 1])
plotmodel(gp)

%% grid search for pre-initialization
fprintf('grid search starts\n');
v = gridsearch(gp,@np_ode_fg);
gp = vec2par(gp,v);
fprintf('Posterior after grid search: %f\n', np_ode_fg(gp));

%% optimization
pars0 = par2vec(gp);
pars = optim_npde(gp,@np_ode_target,50);
gp = vec2par(gp,pars);
fprintf('Posterior after optimization: %f\n', np_ode_fg(gp));
plotmodel(gp)
% wx = 30; wy = wx*0.5;
% set(gcf,'PaperSize',[wx wy], 'PaperPosition',[0 0 wx wy]);
% print('vdp.png','-dpng', '-r300');
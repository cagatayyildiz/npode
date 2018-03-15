% updating the path
addpath(genpath('.'))

% data generation
model = 'vdp'; % Van der Pol oscillator as a toy model
[t,Y] = gen_data('model',model); % time points and noisy data using VDP

% fitting
gp = npode_fit(t,Y); % fits npODE model and returns the learned parameters
% lengthscale and ind. points grid width can be set optionally:
% gp = npode_fit(t,Y,'W',6,'ell',[2 2]); 

% visualization
gp.ode_model = model; % needed to visualize true states
plotmodel(gp) % plots the vector field, true states and trajectories

% predicting future cycles
ts = 0:0.1:40; % new time points
x0 = [2 0]; % initial value
X = npode_predict(gp,ts,x0); % prediction

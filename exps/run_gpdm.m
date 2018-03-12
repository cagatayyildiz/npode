function [X_gpdm,Y_gpdm] = run_gpdm(X,Y,segments,N)
%RUN_GPDM Runs GPDM algorithm.
%   The runner script provided by the authors of the paper is only slightly
%   modified so that correct variables are given as input.
%
%   INPUT
%       X - low dimensional (latent) representation of the training data
%       Y - original training data
%       segments - starting indices of the time series given as input
%       N - number of data points in the data set (training + testing)
%
%   OUTPUT
%       X_gpdm - inferrred latent representation
%       Y_gpdm - reconstructed data

%% check if library exists
global ROOT_GPDM;
if ~exist(ROOT_GPDM,'dir')
    warning('GPDM implementation not found under %s folder.\n',ROOT_GPDM);
    X_gpdm = [];
    Y_gpdm = [];
    return;
end

%% initialization
global USE_GAMMA_PRIOR  % gamma prior for dynamics, only works with RBF kernel
global GAMMA_ALPHA % defines shape of the gamma prior
global USE_LAWRENCE % fix dynamics HPs, as Lawrence suggested (use with thetad = [0.2 0.01 1e6];) 
global FIX_HP % fix all HPs
global MARGINAL_W % marginalize over W while learning X_gpdm
global MARGINAL_DW % marginalize over scale in dynamics while learning X_gpdm
global LEARN_SCALE % use different scales for different output dimensions
global REMOVE_REDUNDANT_SCALE % let W absorb the overall scale of reconstruction
global W_VARIANCE % kappa^2 in the paper, not really the variance though
global M_CONST % M value in Jack's master's thesis
global BALANCE % Constant in front of dynamics term, set to D/q for the B-GPDM
global SUBSET_SIZE % Number of data to select for EM, set -1 for all data. 
global USE_OLD_MISSING_DATA

M_CONST = 1; 
REMOVE_REDUNDANT_SCALE = 1;
LEARN_SCALE = 1; 
MARGINAL_W = 0; 
MARGINAL_DW = 0; 
W_VARIANCE = 1e6; 
FIX_HP = 0; 
USE_GAMMA_PRIOR = 0; 
GAMMA_ALPHA = [5 10 2.5]; 
USE_LAWRENCE = 0;
BALANCE = 1;
SUBSET_SIZE = -1; 

opt = foptions;
opt(1) = 1;
opt(9) = 0;
if MARGINAL_W == 1
    opt(14) = 100; % total number of iterations
    extItr = 1; 
else
    opt(14) = 10; % rescaling every 10 iterations
    extItr = 100; % do extItr*opt(14) iterations in total
end  

% ensure that VGPLVM is not loaded
global ROOT_VGPLVM;
rmpath(genpath(ROOT_VGPLVM));

%% learning
% modelType(1) : input of dynamics
%   0 => [x_t, x_{t-1}]
%   1 => [x_t, x_t - x_{t-1}]
%   2 => [x_t]
% modelType(2) : output of dynamics 
%   0 => x_{t+1} 
%   1 => x_{t+1} - x_t
% modelType(3) : kernel type
%   0 => RBF kernel with weighted dimensions, use with input 0 or 1
%   1 => RBF kernel 
%   2 => Linear kernel
%   3 => weighted Linear kernel + RBF kernel with weighted dimensions, use with
%   input 0 or 1
%   4 => weighted linear kernel
%   5 => linear + RBF

modelType = [2 0 5]; 
theta = [1 1 exp(1)];
thetad = [0.9 1 0.1 exp(1)];
D = size(Y,2);
w = ones(D,1);

[X_gpdm_tr theta thetad w] = gpdmfitFull(X, Y, w, segments, theta, thetad, opt, ... 
     extItr, modelType, []);

%% imputation / forecasting
[K invK] = computeKernel(X_gpdm_tr, theta);
[Xin Xout] = priorIO(X_gpdm_tr, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
simStart = [X_gpdm_tr(segments(1)+1,:) X_gpdm_tr(1,:)]; %  inputs 2 points in case using 2nd order model
[X_gpdm, ~] = simulatedynamics(X_gpdm_tr, segments, thetad, invKd, N, simStart, modelType);

%% reconstruction
KxX = kernel(X_gpdm,X_gpdm_tr, thetaConstrain(theta));
KXX = computeKernel(X_gpdm_tr, theta);
Y_gpdm = KxX / KXX * Y; 
end


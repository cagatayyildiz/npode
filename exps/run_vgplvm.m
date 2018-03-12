function [X_vgplvm,Y_vgplvm] = run_vgplvm(Y,q,N,obs_idx,miss_idx)
%RUN_VGPLVM Runs variational GPLVM algorithm.
%   The runner script provided by the authors of the paper is only slightly
%   modified so that correct variables are given as input.
%
%   INPUT
%       Y - original training data
%       q - number of latent dimensions
%       N - number of data points in the data set (training + testing)
%       obs_idx - indices of the observed data
%       miss_idx - indices of the missing part
%
%   OUTPUT
%       X_vgplvm - inferrred latent representation
%       Y_vgplvm - reconstructed data

%% 
global ROOT_VGPLVM;
if ~exist(ROOT_VGPLVM,'dir')
    warning('VGPLVM implementation not found under %s folder.\n',ROOT_VGPLVM);
    X_vgplvm = [];
    Y_vgplvm = [];
    return;
end
addpath(genpath(ROOT_VGPLVM))

% Define constants
latentDim = q; % this is Q, the number of latent dimensions
indPoints = 100; % number of inducing points
% dynamicKern = {'matern32', 'bias', 'white'};
dynamicKern = {'rbf', 'bias', 'white'}; % kernel k_t for the GP for x(t)
initX ='ppca'; % initialize latent space with ppca

% Corresponding timestamps (artificial and equally spaced for this demo)
t = linspace(0, 2*pi, N+1)';
t = t(1:end-1, 1);
timeStampsTraining = t(obs_idx);
timeStampsTest = t(miss_idx);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'}; % Kernel k_x for the GP for f(x)
options.numActive = indPoints;
options.optimiser = 'scg';

d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
model = vargplvmParamInit(model, model.m, model.X);
model.beta=1/(0.01*var(model.m(:)));
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));

%-------- Add dynamics to the model -----
optionsDyn.type = 'vargpTime';
optionsDyn.t=timeStampsTraining;
optionsDyn.inverseWidth=30;
optionsDyn.initX = initX;

% Dynamic kernel:
kern = kernCreate(t, dynamicKern);
% The following is related to the expected number of
% zero-crossings.(larger inv.width numerator, rougher func)best_gp
if ~strcmp(kern.comp{1}.type,'ou')
    kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
    kern.comp{1}.variance = 1;
end
optionsDyn.kern = kern;

% Fill in with default values whatever is not already set
optionsDyn = vargplvmOptionsDyn(optionsDyn);
model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);

fprintf(1,'# Further calibration of the initial parameters...\n');
model = vargplvmInitDynamics(model,optionsDyn);
model.vardist.parallel=1;
% do not learn beta for few iterations for intitilization
model.learnBeta = 0;
display = 1;
fprintf(1,'# Intitiliazing the model (fixed beta) %d iterations...\n',100);
model = vargplvmOptimise(model, display, 100);
disp('# Saving model after optimising beta...')
% modelWriteResult(model, dataSetName, experimentNo);

%% optimisation
model.learnBeta = 1;
iters = 1000; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
% Save the results.
fprintf(1,'# Saving model after doing %d iterations\n',iters)
% modelWriteResult(model, dataSetName, experimentNo);

% See the final lengthscales (see how some dimensions are switched-off).
% bar(model.kern.comp{1}.inputScales)
 
%% predictions
fprintf('# Only times prediction...\n');
% Prediction using the only information in the test time points
[X_vgplvm, Testcovars] = vargplvmPredictPoint(model.dynamics, t);
Y_vgplvm = vargplvmPosteriorMeanVar(model, X_vgplvm, Testcovars);

rmpath(genpath(ROOT_VGPLVM))
end


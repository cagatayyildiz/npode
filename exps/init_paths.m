function init_paths()
%INIT_PATHS Initializes paths for GPDM, NPODE and VGPLVM
%   You need to update the paths based on your folder structure. Basically,
%   it should be sufficient to update the roots if you download the
%   implementations from the links given below. 
%   Because both GPDM and VGPLVM implementations contain a function named 
%   "computeKernel", the correct library must be loaded. Our code ensures 
%   that by keeping VGPLVM root global and loading/unloading the library on
%   the fly.

%% root directories
global ROOT_VGPLVM;
global ROOT_GPDM;
ROOT_NPODE = '..';
ROOT_GPDM = ''; % update this, see http://www.dgp.toronto.edu/~jmwang/gpdm/
ROOT_VGPLVM = ''; % update this, see https://github.com/SheffieldML/vargplvm

%% paths needed by npode
addpath(fullfile(ROOT_NPODE,'io'));
addpath(fullfile(ROOT_NPODE,'ode'));
addpath(fullfile(ROOT_NPODE,'ode_systems'));
addpath(fullfile(ROOT_NPODE,'plotting'));
addpath(fullfile(ROOT_NPODE,'utils'));

%% paths needed by gpdm
if exist(ROOT_GPDM,'dir')
	addpath(genpath(ROOT_GPDM))
end

end


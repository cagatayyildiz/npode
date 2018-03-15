function demo_cmu_walking(fname,exp_type)
%CMU_WALKING_EXP Runs three algorithms on a CMU walking data sequence
%
% To run the experiments, switch to 'exps' folder.
%
%   Input
%       fname    - name of the amc file to be used in the experiment
%       exp_type - 1 for imputation, 2 for forecasting experiments


if ~exist('fname','var') || isempty(fname)
    fname = '07_02.amc';
end
if ~exist('exp_type','var') || isempty(exp_type)
    exp_type = 1;
end

%% initialize paths for GPDM, NPODE and VGPLVM
% you need to update the paths based on your file structure
init_paths()

%% read data
step = 4; % downsampling rate
q = 3; % number of latent dimensions
[N,Y,X_pca,u,v,X_tr,Y_tr,segments,obs_idx,miss_idx,Y_pca] ...
    = read_data(fname,exp_type,step,q);

%% NPODE
if exp_type == 1
    npde_func = @run_npde_imput;
elseif exp_type == 2
    npde_func = @run_npde_predict;
end
[X_npde,Y_npde] = npde_func(X_tr,u,v,obs_idx,miss_idx);

%% GPDM
[X_gpdm,Y_gpdm] = run_gpdm(X_tr,Y_tr,segments,N);

%% VGPLVM
[X_vgplvm,Y_vgplvm] = run_vgplvm(Y_tr,q,N,obs_idx,miss_idx);

%% plot reconstructions
% plot_reconstr(Y,Y_npde,obs_idx,miss_idx)
% plot_reconstr(Y,Y_gpdm,obs_idx,miss_idx)
% plot_reconstr(Y,Y_vgplvm,obs_idx,miss_idx)

%% errors
pca_obs_err = (Y(obs_idx,:)-Y_pca(obs_idx,:)).^2;
pca_miss_err = (Y(miss_idx,:)-Y_pca(miss_idx,:)).^2;
rmse_pca_obs_err = sqrt(sum(pca_obs_err(:) / numel(pca_obs_err)));
rmse_pca_miss_err = sqrt(sum(pca_miss_err(:) / numel(pca_miss_err)));
fprintf('PCA reconstruction errors for observed and missing parts: %.3f and %.3f\n', ...
    rmse_pca_obs_err, rmse_pca_miss_err);

if exist('Y_npde','var') && numel(Y_npde) == numel(Y)
    np_de_obs_err = (Y(obs_idx,:)-Y_npde(obs_idx,:)).^2;
    np_de_miss_err = (Y(miss_idx,:)-Y_npde(miss_idx,:)).^2;
    rmse_np_de_obs_err = sqrt(sum(np_de_obs_err(:) / numel(np_de_obs_err)));
    rmse_np_de_miss_err = sqrt(sum(np_de_miss_err(:) / numel(np_de_miss_err)));
    fprintf('npDE errors for observed and missing parts: %.3f and %.3f\n', ...
        rmse_np_de_obs_err, rmse_np_de_miss_err);
end

if exist('Y_gpdm','var') && numel(Y_gpdm) == numel(Y)
    gpdm_obs_err = (Y(obs_idx,:)-Y_gpdm(obs_idx,:)).^2;
    gpdm_miss_err = (Y(miss_idx,:)-Y_gpdm(miss_idx,:)).^2;
    rmse_gpdm_obs_err = sqrt(sum(gpdm_obs_err(:) / numel(gpdm_obs_err)));
    rmse_gpdm_miss_err = sqrt(sum(gpdm_miss_err(:) / numel(gpdm_miss_err)));
    fprintf('GPDM errors for observed and missing parts: %.3f and %.3f\n', ...
        rmse_gpdm_obs_err, rmse_gpdm_miss_err);
end

if exist('Y_vgplvm','var') && numel(Y_vgplvm) == numel(Y)
    vgplvm_obs_err = (Y(obs_idx,:)-Y_vgplvm(obs_idx,:)).^2;
    vgplvm_miss_err = (Y(miss_idx,:)-Y_vgplvm(miss_idx,:)).^2;
    rmse_vgplvm_obs_err = sqrt(sum(vgplvm_obs_err(:) / numel(vgplvm_obs_err)));
    rmse_vgplvm_miss_err = sqrt(sum(vgplvm_miss_err(:) / numel(vgplvm_miss_err)));
    fprintf('varGPLVM errors for observed and missing parts: %.3f and %.3f\n', ...
        rmse_vgplvm_obs_err, rmse_vgplvm_miss_err);
end

end






















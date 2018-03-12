function [X_npde,Y_npde]  = run_npde_predict(X,u,v,obs_idx,miss_idx)
%RUN_NPDE_PREDICT Runs npDE algorithm on a data set where the tail is
%missing
%
%   INPUT
%       X - low dimensional (latent) representation of the training data
%       u - principal components
%       v - eigenvalues associated with principal components
%       obs_idx - indices of the observed data
%       miss_idx - indices of the missing part
%
%   OUTPUT
%       X_npode - inferrred latent representation
%       Y_npode - reconstructed data

% introducing artificial time steps
f = 10;
t = obs_idx' / f;
t_all = union(miss_idx,obs_idx)' / f;

% observe that we split the training data into two pieces: the first part,
% denoted by _tr subscript, is the part used for training and the second
% part is set aside for model (lengthscale) selection.
N_tr = floor(length(t) * 3 / 4);
t_tr = t(1:N_tr);
X_tr = X(1:N_tr,:);

%% init handles
f_target = @np_ode_target;
f_post = @np_ode_fg;
f_mean_path = @np_ode_mean_path;

%% create initial points with different length scales
opts.kernel     = 'gauss';
opts.Xlocs      = 'grid';
opts.optpars    = 'log_sn-Fw'; % 'log_sf-log_sn-x0-Fw'
opts.W          = 5;
opts.log_sf     = log(1);
opts.D          = size(X,2);
opts.sn         = (max(X) - min(X)) / 10;
R = 3;
for r = 1:R
    c = floor((r-1)/3)+1;
    if mod(r,1) == 1
        opts.ell    = c*[1 1 1];
    elseif mod(r,1) == 2
        opts.ell    = c*[1 1 0.75];
    else
        opts.ell    = c*[1 0.75 0.75];
    end

    vf_gp = np_de_model({t},{X},opts);
    vf_gp = set_data(vf_gp,{t_tr},{X_tr});
    fprintf('Initial posterior: %f\n',f_post(vf_gp));
    gp0s{r} = vf_gp;
    fprintf('grid search starts\n');
    vec = gridsearch(vf_gp,f_post);
    gp0s{r} = vec2par(vf_gp,vec);
    fprintf('Posterior after grid search: %f\n', f_post(gp0s{r}));
end

%% optimization
M = length(gp0s);
gpts = cell(1,M);
for r = 1:M
    gp_r = gp0s{r};
    gp_r = set_data(gp_r,{t_tr},{X_tr});
    pars = optim_npde(gp_r,f_target,50);
    gp_r = vec2par(gp_r,pars);
    fprintf('Optimized posterior: %f\n', f_post(gp_r));
    gpts{r} = gp_r;
    gp_r
    % plotmodel(gp_r);
end

%% select the model with the best predictive accuracy
best_res = inf;
for i = 1:length(gpts)
    gp_ = gpts{i};
    X_npde = f_mean_path(gp_,{t},gp_.x0);
    X_npde = X_npde{1};
    res = abs(X(N_tr+1:end,:)-X_npde(N_tr+1:end,:)).^2;
    res = sqrt(mean(res(:)));
    if res < best_res
        best_res = res;
        best_gp = gp_;
    end
end
vf_gp = best_gp;
fprintf('Selected posterior: %f\n',f_post(vf_gp));
% plotmodel(vf_gp);

%% re-train the best model
vf_gp = set_data(vf_gp,{t},{X});
pars = optim_npde(vf_gp,f_target,100);
vf_gp = vec2par(vf_gp,pars);
fprintf('posterior at iter %d is %f\n',c,f_post(vf_gp));

%% reconstruction
X_npde = f_mean_path(vf_gp,{t_all},vf_gp.x0);
X_npde = X_npde{1};
q = size(X,2);
Y_npde = X_npde * diag(sqrt(v(1:q))) * u(:, 1:q)';

end






















function [X_npode,Y_npode] = run_npde_imput(X,u,v,obs_idx,miss_idx)
%RUN_NPDE_IMPUT Runs npODE algorithm on a data set where the missing data
%is in the middle.
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

ts{1} = t(1:miss_idx(1)-1);
Xs{1} = X(1:miss_idx(1)-1,:);
ts{2} = t(miss_idx(1):end);
Xs{2} = X(miss_idx(1):end,:);

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

    vf_gp = np_de_model({t(1:miss_idx(1)-1)},{X(1:miss_idx(1)-1,:)},opts);
    
    fprintf('Initial posterior: %f\n',f_post(vf_gp));
    pars = par2vec(vf_gp);

    pars = gridsearch(vf_gp,f_post);
    gp0s{r} = vec2par(vf_gp,pars);
    fprintf('Posterior after grid search: %f\n', f_post(gp0s{r}));

    pars = optim_npde(vf_gp,f_target,50,pars);
    gp0s{r} = vec2par(vf_gp,pars);
    fprintf('Posterior after optimization: %f\n', f_post(gp0s{r}));
end

%% sort initial points based on predictive accuracy
preds = zeros(1,length(gp0s));
for i = 1:length(gp0s)
    X_npode = f_mean_path(gp0s{i},{t_all(obs_idx)},gp0s{i}.x0);
    X_npode = X_npode{1};
    seg2_start = miss_idx(end)+1;
    res = abs(X(seg2_start:end,:)-X_npode(seg2_start:end,:)).^2;
    preds(i) = sqrt(mean(res(:)));
end
[~,idx] = sort(preds);
gp0s = gp0s(idx);

%% training: from the scratch
% vf_gp = gp0s{1};
% X_filled = f_mean_path(vf_gp,{t_all},vf_gp.x0);
% X_filled = X_filled{1};
% X_filled(obs_idx,:) = X;

% use the first gp
vf_gp = gp0s{1};
vf_gp = set_data(vf_gp,{t},{X});

pars = optim_npde(vf_gp,f_target,100);
vf_gp = vec2par(vf_gp,pars);
% plotmodel(vf_gp)

% training error
X_npode = f_mean_path(vf_gp,{t_all(obs_idx)},vf_gp.x0);
X_npode = X_npode{1};
% res = abs(X-X_npode).^2;
% sqrt(mean(res(:)))


%% reconstruction
X_npode = f_mean_path(vf_gp,{t_all},vf_gp.x0);
X_npode = X_npode{1};
q = size(X,2);
Y_npode = X_npode * diag(sqrt(v(1:q))) * u(:, 1:q)';

end






















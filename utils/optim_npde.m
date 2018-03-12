function [pars,vf_gp] = optim_npde(vf_gp,ftarget,num_iter,pars,batch)
%OPTIM_NP_DE Optimizes npODE model using L-BFGS / fminunc
%
%   INPUT
%       vf_gp - np_DE model to be optimized
%       ftarget - target function to be minimized
%       num_iter - optional num. of iter. of optimization (default 250 for
%                  batch opt. and 25 for incremental opt)
%       pars - optional initial point vector (default is vf_gp parameters)
%       batch - optional boolean indicating batch/incremental optimization
%               (default true)
%
%   OUTPUT
%       pars - optimized parameters

%% init the params
if ~exist('pars','var') || isempty(pars)
    pars = par2vec(vf_gp);
end  
if ~exist('batch','var')
    batch = true;
end  
if ~exist('num_iter','var') || isempty(num_iter)
    if batch
        num_iter = 250;
    else
        num_iter = 25;
    end
end  

%% init optimizer and its params
lbfgs = false;
if lbfgs
    opt_func = @minFunc;
    opt_args.maxIter = num_iter;
    opt_args.maxFunEvals = num_iter;
    opt_args.display = 'full';
    opt_args.progTol = 1e-12;
    opt_args.DerivativeCheck =  0;
    opt_args.Method = 'lbfgs'; % 'lbfgs' 'sd' 'csd' 'bb' 'newton0' 'pnewton0' 'cg' 'scg' 'pcg'
else
    opt_func = @fminunc;
    opt_args = optimoptions('fminunc','Display','iter', ...
        'OptimalityTolerance',1e-12, ...
        'SpecifyObjectiveGradient',true, ...
        'MaxIterations',num_iter);
end

%% run the optimization
if batch % batch mode
    try
        pars = opt_func(@(par)ftarget(vf_gp,par), pars, opt_args);
        vf_gp = vec2par(vf_gp,pars);
    catch ex
        fprintf('Optimization stopped with an error: %s\n',ex.message);
    end
else % incremental optimization
    t = vf_gp.t;
    Y = vf_gp.Y;
    step = 2;
    C = length(t{1})/step; % number of batches 
    for c = 1:C
        sz = step*c;
        for i = 1:vf_gp.Nt
            Y_{i} = Y{i}(1:sz,:);
            t_{i} = t{i}(1:sz);
        end
        vf_gp = set_data(vf_gp,t_,Y_);
        try
            pars = opt_func(@(par)ftarget(vf_gp,par), pars, opt_args);
            vf_gp = vec2par(vf_gp,pars);
        catch ex
            fprintf('Optimization stopped with an error: %s\n',ex.message);
        end
        fprintf('value of the target fnc at batch %d is %f\n',c,ftarget(vf_gp,pars));
    end
end

end


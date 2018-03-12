
function [D,x0,t,X,Y,ode_fun] = gen_data(opts)

%% set defaults
if ~isfield(opts,'Ny')
    opts.Ny = 25;
end
if ~isfield(opts,'t_ends')
    if strcmp(opts.model,'vdp')
        opts.t_ends = 7;
    elseif strcmp(opts.model,'fhn')
        opts.t_ends = 10;
    elseif strcmp(opts.model,'lv')
        opts.t_ends = 5;
    elseif strcmp(opts.model,'calcium')
        opts.t_ends = 7;
    end
end
if ~isfield(opts,'t_starts')
    opts.t_starts = 0.01;
end
if ~isfield(opts,'draw')
    opts.draw = false;
end
if ~isfield(opts,'t_slice')
    opts.t_slice = 'uniform';
end

if isscalar(opts.Ny)
    opts.Ny = opts.Ny * ones(opts.Nt,1);
end
if isscalar(opts.t_ends)
    opts.t_ends = opts.t_ends * ones(opts.Nt,1);
end
if isscalar(opts.t_ends)
    opts.t_ends = opts.t_ends * ones(opts.Nt,1);
end
if isscalar(opts.t_starts)
    opts.t_starts = opts.t_starts * ones(opts.Nt,1);
end

%% set ode function
if strcmp(opts.model,'vdp')
    ode_fun = @(t,x) ode_vdp(t,x);
    D = 2;
elseif strcmp(opts.model,'fhn')
    ode_fun = @(t,x) ode_fhn(t,x);
    D = 2;
elseif strcmp(opts.model,'lv')
    ode_fun = @(t,x) ode_lv(t,x);
    D = 2;
elseif strcmp(opts.model,'calcium')
    ode_fun = @(t,x) ode_calcium(t,x);
    D = 4;
end


%% generate data
X = cell(opts.Nt,1);
Y = cell(opts.Nt,1);
t = cell(opts.Nt,1);

% initial points
if ~isfield(opts,'x0')
    tmp = lhsdesign(opts.Nt,D);
    % plot(tmp(:,1), tmp(:,2), '*')
    if strcmp(opts.model,'vdp')
        x0 = tmp .* [4,4] - [2,2];
    elseif strcmp(opts.model,'fhn')
%         x0 = tmp .* [5,2.5] - [2,1];
        x0 = tmp .* [6,4] - [3,2];
    elseif strcmp(opts.model,'lv')
        x0 = tmp .* [4,4] + [6,1];
    elseif strcmp(opts.model,'calcium')
        x0 = tmp .* [4,4] - [2,2];
    end
else
    x0 = opts.x0;
end

% data itself
for i = 1:opts.Nt
    if strcmp(opts.t_slice,'uniform')
        t{i} = linspace(opts.t_starts(i),opts.t_ends(i),opts.Ny(i))';
    elseif strcmp(opts.t_slice,'random')
        tmp = cumsum(rand(opts.Ny(i),1));
        t{i} = opts.t_ends(i) * tmp ./ tmp(end,:);
    end

    [~,X_i] = ode45(ode_fun, t{i}, x0(i,:));
    X{i} = X_i;
    Y{i} = X_i + opts.sn_true*randn(size(X_i));
end

%% visualize
if opts.draw
    figure;
    for i = 1:data_opts.Nt
        subplot(ceil(sqrt(opts.Nt)),ceil(sqrt(opts.Nt),i))
        plot(t{i},Y{i}(:,1))
        hold on
        plot(t{i},Y{i}(:,2))
        hold off
    end
end

end
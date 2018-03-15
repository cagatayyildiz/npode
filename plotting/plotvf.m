function [] = plotvf(gp,idx)
% PLOTVF Plots the vector field and trajectories
%
% INPUT
%       gp  - npODE instance storing underlying parameters
%       idx - the index of the input time series to be plotted

if size(gp.X,2) ~= 2
    return;
end

if ~exist('idx','var')
    idx = 1:gp.Nt;
end

qscale = 0.1;
plot(nan);
hold on;

%% vector field
% on a grid
Ng = 20;
d1 = (max(gp.X(:,1)) - min(gp.X(:,1))) / 6;
d2 = (max(gp.X(:,2)) - min(gp.X(:,2))) / 6;
v1 = linspace(min(gp.X(:,1))-d1, max(gp.X(:,1))+d1, Ng);
v2 = linspace(min(gp.X(:,2))-d2, max(gp.X(:,2))+d2, Ng);
[X1,X2] = meshgrid(v1, v2);
Xs = [X1(:) X2(:)];
Fg = predf_scalar(gp, Xs);

h = quiver(X1(:),X2(:),qscale*Fg(:,1),qscale*Fg(:,2), 'autoscale','off','color','k', 'linewidth',1.0);
set(h, 'color',0.7*[1 1 1]);
% quiverc(X1(:),X2(:),qscale*Fg(:,1),qscale*Fg(:,2));
ip = plot(gp.X(:,1),gp.X(:,2),'ko', 'linewidth',1.5);

%% inducing vectors
iv = quiver(gp.X(:,1),gp.X(:,2),qscale*gp.F(:,1),qscale*gp.F(:,2),'autoscale','off','color','k', 'linewidth',1.5);

%% data points
for i = idx
%     dp = plot(gp.Y{i}(:,1), gp.Y{i}(:,2), strcat('bo'), 'markersize',7, 'linewidth',2.0);
    dp = plot(gp.Y{i}(:,1), gp.Y{i}(:,2), strcat(gp.clrs(i),'o'), 'markersize',7, 'linewidth',2.0);
%     plot(gp.Y{i}(:,1), gp.Y{i}(:,2), strcat(gp.clrs(i),':'), 'markersize',5, 'linewidth',1.0);
end

%% trajectories + initial values
ts = cell(1,length(idx));
for i = idx
    ts{i} = linspace(0.001,gp.t{i}(end),200)';
end
Xp = np_ode_mean_path(gp,ts,gp.x0);
for i = idx
    tr = plot(Xp{i}(:,1), Xp{i}(:,2), strcat(gp.clrs(i),'-'), 'linewidth', 1.5);
    plot(gp.x0(i,1),gp.x0(i,2),'o','markersize',5,'linewidth',5,'Color', gp.clrs(i))
%     plot(Xp{i}(:,1), Xp{i}(:,2), strcat(gp.clrs(i),'-'), 'markersize',5, 'linewidth',1.0);
    if ~isempty(gp.x0_true) 
        x0_true = gp.x0_true(i,:);
    else
        x0_true = gp.Y{i}(1,:);
    end
    if ~isempty(gp.x0_true) && ~isempty(gp.ode_fun)
        [~,Xsol] = ode45(gp.ode_fun,ts{i},x0_true);
        plot(gp.x0_true(i,1),gp.x0_true(i,2),'*','linewidth',5,'Color', gp.clrs(i))
        true_tr = plot(Xsol(:,1), Xsol(:,2),'--', 'linewidth',1.5, 'Color', gp.clrs(i));
    end
end


%% finish
hold off;
xlim([min(v1) max(v1)]);
ylim([min(v2) max(v2)]);
xlabel('$x_1$', 'interpreter','latex');
ylabel('$x_2$', 'interpreter','latex');
if exist('true_tr','var')
    legend([dp,tr,true_tr,ip,iv], ...
    {"data points","estim. traj","true traj","ind. points","ind. vectors"})
else
    legend([dp,tr,ip,iv], ...
    {"data points","estim. traj","ind. points","ind. vectors"})
end
title(sprintf('[logp %.2f] [sf %.2f] [ell %.3f %.3f] [sn %.3f %.3f]', ...
    np_ode_fg(gp),gp.sf,gp.ell(1), gp.ell(2), gp.sn(1), gp.sn(2)));


end

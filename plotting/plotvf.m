% plot the ODE system vector field
%
% only for two-dimensional problems
%
% truef: plot also true function
%
function [] = plotvf(gp,x0,truef)

if ~exist('x0','var')
    x0 = [];
end
if ~exist('truef','var')
    truef = [];
end

if size(gp.X,2) ~= 2
    return;
end

qscale = 0.1;
plot(nan);
hold on;

%% data points
for i = 1:gp.Nt
    dp = plot(gp.Y{i}(:,1), gp.Y{i}(:,2), strcat(gp.clrs(i),'o'), 'markersize',7, 'linewidth',2.0);
%     plot(gp.Y{i}(:,1), gp.Y{i}(:,2), strcat(gp.clrs(i),':'), 'markersize',5, 'linewidth',1.0);
end

%% inducing vectors
iv = quiver(gp.X(:,1),gp.X(:,2),qscale*gp.F(:,1),qscale*gp.F(:,2),'autoscale','off','color','k', 'linewidth',1.5);

%% trajectories + initial values
ts = cell(1,gp.Nt);
for i = 1:gp.Nt
    ts{i} = linspace(0.001,gp.t{i}(end),200)';
end
Xp = np_ode_sens(ts,gp.x0,gp.F,gp.X,gp.ell,gp.sf);
for i = 1:gp.Nt
    tr = plot(Xp{i}(:,1), Xp{i}(:,2), strcat(gp.clrs(i),'-'), 'linewidth', 1.5);
    plot(gp.x0(i,1),gp.x0(i,2),'o','markersize',5,'linewidth',5,'Color', gp.clrs(i))
%     plot(Xp{i}(:,1), Xp{i}(:,2), strcat(gp.clrs(i),'-'), 'markersize',5, 'linewidth',1.0);
    if ~isempty(gp.x0_true)
        [~,Xsol] = ode45(gp.ode_fun,ts{i},gp.x0_true(i,:));
        plot(gp.x0_true(i,1),gp.x0_true(i,2),'*','linewidth',5,'Color', gp.clrs(i))
        true_tr = plot(Xsol(:,1), Xsol(:,2),'--', 'linewidth',1.5, 'Color', gp.clrs(i));
    end
end

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


%% error lines
for i = 1:gp.Nt
%     line([Xp{i}(:,1),gp.Y{i}(:,1)]',[Xp{i}(:,2),gp.Y{i}(:,2)]', 'Color', gp.clrs(i))
end

%% true vf
if ~isempty(truef)
    [~,truex] = ode45(truef,ts,x0);
    plot(truex(:,1), truex(:,2), 'r--', 'linewidth',2);
end

%% finish
hold off;
xlim([min(v1) max(v1)]);
ylim([min(v2) max(v2)]);
xlabel('Var 1');
ylabel('Var 2');
legend([dp,tr,true_tr,ip,iv],{"data points","est. traj","true traj","ind. points","ind. vectors"})
title(sprintf('[logp %.2f] [sf %.2f] [ell %.3f %.3f] [sn %.3f %.3f]',np_ode_fg(gp),gp.sf,gp.ell(1), gp.ell(2), gp.sn(1), gp.sn(2)));

end

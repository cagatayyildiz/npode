function [] = plotode(gp,s)
% PLOTODE Plots the trajectories (true and inferred) and uncertainty
%
% INPUT
%       gp  - npODE instance storing underlying parameters
%       s   - the index of the input time series to be plotted


% grid points
D = gp.D;
colors = ggplotcolors(D);

plot(nan);
hold on;

for i=1:D
    dp(i) = plot(gp.t{s}, gp.Y{s}(:,i), '.', 'markersize',15, 'color', colors(i,:));
end

if D == 2
    Ns = 600;
    ts = linspace(0.001,4*max(gp.t{s}),Ns)';
else
    Ns = 100;
    ts = linspace(0.001,max(gp.t{s}),Ns)';
end
Xnpde = np_ode_mean_path(gp,{ts},gp.x0(s,:));
Xnpde = Xnpde{1};
for i=1:D
    sig(i) = fill([ts; flip(ts)], [Xnpde(:,i)-2*gp.sn(i); flip(Xnpde(:,i)+2*gp.sn(i))], ...
        colors(i,:),'facealpha',0.40,'edgecolor','none');
end
if ~isempty(gp.x0_true) 
    x0 = gp.x0_true(s,:);
else
    x0 = gp.Y{s}(1,:);
end
if ~isempty(gp.ode_fun)
    [~,Xmu] = ode45(gp.ode_fun,ts,x0);
    for i=1:D
        mp(i) = plot(ts,Xmu(:,i),'k-');
    end
end

hold off;
ylim([min(gp.Y{s}(:))-1 max(gp.Y{s}(:))+1]);
xlim([0 max(ts)]);
xlabel('time');
ylabel('x(t)');
title(sprintf('States'))
if D == 2
    if s == 1
        if exist('mp','var')
        legend([mp(1),sig,dp],{"True $x(t)$","Estimated $x_1(t)+2\sigma_n$",...
            "Estimated $x_2(t)+2\sigma_n$","Data $y_1$","Data $y_2$"},...
            'interpreter','latex');
        else
        legend([sig,dp],{"Estimated $x_1(t)+2\sigma_n$",...
            "Estimated $x_2(t)+2\sigma_n$","Data $y_1$","Data $y_2$"},...
            'interpreter','latex');
        end
    end
end
    
end



% plot the ODE system
%
% truef: plot also true function
%
function [] = plotode(gp,s,x0,truef,Nsamples)
	
	if ~exist('truef','var')
		truef = [];
	end
	if ~exist('Nsamples','var')
		Nsamples = 0;
	end

	% grid points
	Ns = 200;
    
    t_ = [];
    for i=1:gp.Nt, t_ = [t_;gp.t{i}]; end
	ts = linspace(0.001,max(t_),Ns)';
	D = gp.D;
	
	colors = ggplotcolors(D);

	% integrate the system
% 	[Xmu,Xcov,Xz] = num_int1(gp,x0,ts,Nsamples);
% 	Xsd = zeros(Ns,D);
% 	Xysd = zeros(Ns,D);
% 	for i=1:D
% 		Xsd(:,i) = sqrt(squeeze(Xcov(i,i,:)));
% 		Xysd(:,i) = sqrt(squeeze(Xcov(i,i,:)) + gp.sn(i)^2);
% 	end
	
    Xmu = np_ode_sens({ts},gp.Y{s}(1,:),gp.F,gp.X,gp.ell,gp.sf,gp.tol);
    Xmu = Xmu{1};
    
	plot(nan);
	hold on;

	for i=1:D
		plot(gp.t{s}, gp.Y{s}(:,i), '.', 'markersize',15, 'color', colors(i,:));
	end
	
	if ~isempty(truef)
		[ts,truex] = ode45(truef,ts,x0);
		for i=1:D
			plot(ts,truex(:,i), '-', 'color', colors(i,:), 'linewidth',2);
		end
    end
    
	if Nsamples > 0
		for i=1:D
			plot(ts, squeeze(Xz(:,i,:)), '-', 'color', [1 1 1]*0.5, 'linewidth',0.5);
		end
	end
	
%	for i=1:D
%		plot(ts, Xmu);
%	end

	for i=1:D
%		fill([ts; flip(ts)], [Xmu(:,i)-2*Xsd(:,i); flip(Xmu(:,i)+2*Xsd(:,i))], colors(i,:),'facealpha',0.30,'edgecolor','none');
%		fill([ts; flip(ts)], [Xmu(:,i)-2*Xysd(:,i); flip(Xmu(:,i)+2*Xysd(:,i))], colors(i,:),'facealpha',0.30,'edgecolor','none');
		fill([ts; flip(ts)], [Xmu(:,i)-2*gp.sn(i); flip(Xmu(:,i)+2*gp.sn(i))], colors(i,:),'facealpha',0.40,'edgecolor','none');
	end
	
	% plot force directions
%	for i=1:D
%		I = gp.F(:,i)>0;
%		tp = gp.t(I);
%		tn = gp.t(~I);
%		xip = gp.X(I,i);
%		xin = gp.X(~I,i);
%		fip = (gp.F(I,i))*(1/5);
%		fin = -abs(gp.F(~I,i))*(1/5);
%		plot([tp tp]', [xip xip+fip]', 'k-', tp, xip+fip, 'k^');%, tp, xip, 'ko');
%		plot([tn tn]', [xin xin+fin]', 'k-', tn, xin+fin, 'kv');%, tn, xin, 'ko');
%	end
	
	hold off;
	ylim([min(gp.Y{s}(:))-1 max(gp.Y{s}(:))+1]);
	xlim([0 max(gp.t{s})]);
	xlabel('time');
	ylabel('x(t)');
    title(sprintf('Trajectory-%d',s))

% 	title(sprintf('[sf %.2f] [ell %.3f %.3f] [sn %.3f %.3f]',gp.sf,gp.ell(1), gp.ell(2), gp.sn(1), gp.sn(2)));
end



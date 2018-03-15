function [f,g,Xmu] = np_ode_fg(gp,par,Xmu)
% NP_ODE_FG Computes the value and the gradients of the ODE posterior
%
% INPUT
%       gp 
%       pars - optional parameters

% OUTPUT
%       f - log-posterior
%       g - gradients of the log posterior
%       Xmu - mean path

if ~exist('par','var')
    par = [];
end
if ~isempty(par)
    gp = vec2par(gp,par);
end
if ~exist('Xmu','var') 
    if nargout == 1 % just the path
        Xmu = np_ode_sens(gp.t,gp.x0,gp.F,gp.X,gp.ell,gp.sf);
    elseif nargout > 1 % path + sensitivities
        [Xmu, S] = np_ode_sens(gp.t,gp.x0,gp.F,gp.X,gp.ell,gp.sf);
    end
end
% gp

%% log-posterior
fth = log_gam(gp.sf,gp.sf_alpha,gp.sf_beta) + ...
      sum(log_gam(gp.sn,gp.sn_alpha,gp.sn_beta)) + ...
      sum(log_gam(gp.ell,gp.ell_alpha,gp.ell_beta)) ;
ff = sum(diag(logmvnpdf(gp.F, 0, gp.Kf)));
fy = 0;
for i = 1:gp.Nt
    lls = cellfun( @(s) logmvnpdf(gp.Y{i}(s,:)', Xmu{i}(s,:)', diag(gp.sn.^2)), num2cell(1:gp.Ny(i))');
    fy = fy + sum(lls);
end
f = ff + fy + fth;

%% gradients
if nargout > 1
    
    dx0 = []; dFw = []; dlog_sf = []; dlog_sn = []; dlog_ell = [];
    
%     non_cent_ = gp.non_cent;
%     gp.non_cent = true;
    
    if contains(gp.optpars,'log_sf')
        dlog_sf = grad_logsf(gp);
    end
    
    if contains(gp.optpars,'log_sn')
        dlog_sn = grad_logsn(gp,Xmu);
    end
    
    if contains(gp.optpars,'log_ell')
        dlog_ell = grad_logell(gp);
    end
    
    if contains(gp.optpars,'x0')
        dx0 = grad_x0(gp,Xmu,S);
    end
    
    if ismember('Fw',gp.optpars)
        dFw = grad_Fw(gp,Xmu,S);
    end
    
%     gp.non_cent = non_cent_;
    
    g = wrap_params(gp,'x0',dx0,'log_sf',dlog_sf,'log_sn',dlog_sn,...
        'log_ell',dlog_ell,'Fw',dFw);
    
    %     plotmodel(gp);
    %     drawnow;
end
end


%% grad log_sn
function dlog_sn = grad_logsn( gp,Xmu )
if ~exist('Xmu','var') 
    [Xmu, ~] = np_ode_sens(gp.t,gp.x0,gp.F,gp.X,gp.ell,gp.sf);
end
dsn_lhood = zeros(size(gp.sn));
for i = 1:gp.Nt
    err = (gp.Y{i}-Xmu{i});
    dsn_lhood = dsn_lhood + sum(err.^2,1)./gp.sn.^3 - gp.Ny(i)./gp.sn;
end
[~,dsn_pri] = log_gam(gp.sn,gp.sn_alpha,gp.sn_beta);
dsn = dsn_lhood + dsn_pri;
dlog_sn = dsn .* gp.sn;
end


%% grad log_sf
function dlog_sf = grad_logsf( gp )
delta = 1e-5;
vec = par2vec(gp);

vec(1) = vec(1) + delta; 
f_sfp = np_ode_fg(gp,vec);

vec(1) = vec(1) - 2*delta; 
f_sfn = np_ode_fg(gp,vec);

dlog_sf = (f_sfp-f_sfn) / (2*delta);
end


%% grad log_ell
function dlog_ell = grad_logell( gp )
delta = 1e-5;

i = 0;
if contains(gp.optpars,'log_sf')
    i = i+1;
end
if contains(gp.optpars,'log_sn')
    i = i+gp.D;
end

vec = par2vec(gp);
dlog_ell = zeros(1,gp.D);
for d = 1:gp.D
    vec(i+d) = vec(i+d) + delta; 
    fp = np_ode_fg(gp,vec);
    vec(i+d) = vec(i+d) - 2*delta; 
    fn = np_ode_fg(gp,vec);
    dlog_ell(d) = (fp-fn) / (2*delta);
    vec(i+d) = vec(i+d) + delta; 
end

end


%% grad Fw
function dFw = grad_Fw( gp,Xmu,S )
if ~exist('Xmu','var') 
    [Xmu, S] = np_ode_sens(gp.t,gp.x0,gp.F,gp.X,gp.ell,gp.sf);
end
dFf = -gp.Kf\gp.F;
dFy = zeros(size(dFf));
for i = 1:gp.Nt
    dFy_i = grad_sens(gp,gp.Y{i},Xmu{i},S{i});
    dFy = dFy + reshape(sum(dFy_i,2),gp.Nf,gp.D); 
end
dF = dFf + dFy;
dFw = gp.L' * dF; % whiten
end

% computation of grad F using sensitivities
function dFy = grad_sens(gp,Y,Xmu,S)
Ny = size(Y,1);
dFy = zeros(gp.Nf*gp.D,Ny);
for i = 1:gp.Nf*gp.D % derivative wrt \theta_i
    for t_ = 1:Ny
        dFy(i,t_) = sum( (Y(t_,:)-Xmu(t_,:)) .* S(:,i,t_)' ./ gp.sn.^2 );
    end
end
end


%% grad x0
function dx0 = grad_x0( gp,Xmu,S )

if ~exist('Xmu','var') 
    [Xmu, S] = np_ode_sens(gp.t,gp.x0,gp.F,gp.X,gp.ell,gp.sf);
end

dx0 = zeros(gp.Nt,gp.D);

for i = 1:gp.Nt % derivative wrt i'th initial value: x_0(i,:)
    err = (gp.Y{i}-Xmu{i}) ./ gp.sn.^2;    
    for j = 1:gp.D % derivative wrt j'th dimension:  x_0(i,j)
        partial_der = squeeze(S{i}(:,numel(gp.F)+j,:))';
        grad = err .* partial_der;
        dx0(i,j) = sum( grad(:) );
    end
end

end
function [X, S] = np_ode_sens(ts,x0,U,Z,ell,sf,tol)
% NP_ODE_SENS Computes the mean path as well as sentitivities
%
% INPUT
%       ts  - time indices (cell array)
%       x0  - initial values (rows of x0)
%       U   - inducing vectors
%       Z   - inducing point locations
%       ell - lengthscale
%       sf  - signal variance
%       tol - tolerance (optional)
%
% OUTPUT
%       X   - mean path
%       S   - sensitivities

if ~exist('tol','var') || isempty(tol)
    tol = 0.01;
end
Nt = length(ts);
X = cell(Nt,1);
S = cell(Nt,1);

for i = 1:Nt
    t_i = ts{i};
    x0_i = x0(i,:);
    if nargout == 1
        X_i = ode_sens(x0_i, t_i, U, Z, ell, sf, tol);
    elseif nargout == 2
        [X_i,S_i] = ode_sens(x0_i, t_i, U, Z, ell, sf, tol);
        S{i} = S_i;
    end
    X{i} = X_i;
end

end

function [Xs,dU] = ode_sens(x0, ts, U, Z, ell, sf, sz)
% ODE_SENS Computes the mean path as well as sentitivities for a single
% time series input
%
% INPUT
%       x0  - initial values (row vector)
%       ts  - time indices (vector)
%       U   - inducing vectors
%       Z   - inducing point locations
%       ell - lengthscale
%       sf  - signal variance
%       sz  - tolerance (optional)
%
% OUTPUT
%       Xs  - mean path
%       dU  - sensitivities

if ~exist('sz','var') || isempty(sz)
    sz = 0.01;
end

Nt = length(ts);
[N,D] = size(x0);
M = size(Z,1);
MD = D*M;
zMDD = [zeros(D,MD) eye(D)];
% 	zMDD = zeros(D,MD);

if nargout == 1
    [~,Xs] = ode45( @(t,x) ode_kernel(x', U, Z, ell, sf, sz)', ts, x0');
else
    aug0 = [x0(:); zMDD(:)];
    [~,augs] = ode45( @(t,aug) ode_kernelsens(aug, U, Z, ell, sf, sz), ts, aug0);

    Xs = augs(:,1:D);
    dU = augs(:,D+1:end);
%   dU = reshape(dU', D, MD, Nt);
    dU = reshape(dU', D, MD+D, Nt);
end
end


function dX = ode_kernel(X, U, Z, ell, sf, sz)
% ODE_KERNEL ODE differential function 
% used only when the mean path is computed, doesnt deal with sensitivities
%
% INPUT
%       X   - states (in rows)
%       U   - inducing vectors
%       Z   - inducing point locations
%       ell - lengthscale
%       sf  - signal variance
%       sz  - tolerance
%
% OUTPUT
%       dX  - differential function values
KxZ = gausskernel(X,Z,ell,sf);
KZZ = gausskernel(Z,Z,ell,sf,sz);
dX = KxZ/KZZ*U;
end


function dxs = ode_kernelsens(state, U, Z, ell, sf, sz)
% ODE_KERNELSENS ODE differential function 
% used whenever sensitivities are computed
%
% INPUT
%       state   - state + sensitivities
%       U       - inducing vectors
%       Z       - inducing point locations
%       ell     - lengthscale
%       sf      - signal variance
%       sz      - tolerance
%
% OUTPUT
%       dxs     - differential computed at state
	
[M,D] = size(U);
MD = M*D;

x = state(1:D);
s = state(D+1:end);
% 	S = reshape(s,D,MD); % MD for U, D for x0
S = reshape(s,D,(MD+D)); % MD for U, D for x0

% state diff
KxZ = gausskernel(x',Z,ell,sf);
KZZ = gausskernel(Z,Z,ell,sf,sz);
invK = inv(KZZ);
dX = KxZ*invK*U;
dx = dX(:);

% 	J = (- 1/ell.^2 * ((x - Z') .* KxZ) * invK * U)';  
J = (- ((x - Z') .* KxZ)./(ell.^2)' * invK * U)';  % vectorised form
R =  [kron(eye(D), KxZ*invK) zeros(D)];            % vectorised form
% 	R =  [kron(eye(D), KxZ*invK)];            % vectorised form

dS = J * S + R;
ds = dS(:);

dxs = [dx; ds];
end





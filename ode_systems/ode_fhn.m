% Fitz-Hugh-Nagumo ODE model
% 
% to solve:
%  [ts,X] = ode45(@(t,x) ode_fhn(x),ts,x0);
%
function dx = ode_fhn(t,x,p)
	if ~exist('p','var') % default parameter values
		p = [0.2,0.2,3]';
	end
	if ~exist('x','var') % default initial values
		x = [-1 1]';
	end

	dx = [p(3).*(x(1,:)-x(1,:).^3./3+x(2,:));
			-1./p(3).*(x(1,:)-p(1)+p(2)*x(2,:))];
end



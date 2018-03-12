% The Calcium ODE model
%
% states (4): Galpha, PLC, CaCyt, CaEr
% parameters (17): k1...k17
% init values: [0.12 0.31 0.0058 4.3]
% equations:
%   Galpha' ==  k(1) + k(2)*Galpha - k(3)*PLC*f(Galpha,K(12)) -  k(4)*CaCyt*f(Galpha,K(13));
%   PLC'    ==  k(5) * Galpha - k(6)*f(PLC,K(14));
%   CaCyt'  ==  k(7) * PLC * CaCyt * f(CaEr,K(15)) + k(8)*PLC+k(9)*Galpha - k(10)*f(CaCyt,K(16)) - k(11)*f(CaCyt,K(17));
%   CaEr'   == -k(7) * PLC * CaCyt * f(CaEr,K(15)) + k(11)*f(CaCyt,K(17));
%
% to solve:
%  [ts,X] = ode45(@(t,x) ode_calcium(x),ts,x0);
%
function dx = ode_calcium(t,x,p) 
	if ~exist('p','var')
		p = [0.09 2 1.27 3.73 1.27 32.24 2 0.05 13.58 153 4.85 0.19 0.73 29.09 2.67 0.16 0.05]';
	end
	if ~exist('x','var')
		x = [0.12 0.31 0.0058 4.3]';
	end

	f = @(x,p) x./(x+p);

	dx = [p(1) + p(2).*x(1,:) - p(3).*x(2,:).*f(x(1,:),p(12)) - p(4).*x(3,:).*f(x(1,:),p(13));
			p(5) .* x(1,:) - p(6).*f(x(2,:),p(14));
			p(7) .* x(2,:) .* x(3,:) .* f(x(4,:),p(15)) + p(8).*x(2,:)+p(9).*x(1,:) - p(10).*f(x(3,:),p(16)) - p(11).*f(x(3,:),p(17));
			-p(7) .* x(2,:) .* x(3,:) .* f(x(4,:),p(15)) + p(11).*f(x(3,:),p(17))];
end



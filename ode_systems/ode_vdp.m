% van der Pol ODE system
%
% to run:
%  [ts,V] = ode45(@(t,y) vdp(t,y),ts,[-1 1]);
%
function dx = ode_vdp(t,x)
	
	% default initial value
	if ~exist('x','var')
		x = [2 0];
	end

	dx = [x(2); 
		 (1-x(1).^2).*x(2)-x(1)];
%      dx = [x(:,2) (1-x(:,1).^2).*x(:,2)-x(:,1)];
end



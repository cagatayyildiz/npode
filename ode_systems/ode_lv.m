% Lotka-Volterra ODE
%
function dx = ode_lv(t,x,p)
	
	if ~exist('p','var')
		p = [2 1 4 1];
	end

	dx = [p(1)*x(1) - p(2)*x(1)*x(2);
		  -p(3)*x(2) + p(4)*x(1)*x(2)];
end


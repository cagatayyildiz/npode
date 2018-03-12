%
% interpolate GP models
%
function [Fmu,Fsd,Fcov] = preddg(gp, ts)

	Ns = length(ts);

	Fmu = zeros(Ns,gp.D);
	Fsd = zeros(Ns,gp.D);
	Fcov = cell(gp.D,1);
	
	for i=1:gp.D
		y = gp.Y(:,i);
		
		Ky = gausskernel(gp.t,gp.t,gp.ell(i),gp.sf(i),gp.sn(i));

		% derivatives
		Dsx = gausskernel(ts,gp.t,gp.ell(i),gp.sf(i),0,1);
		Dss = gausskernel(ts,ts,gp.ell(i),gp.sf(i),0,2);

		dA = Dsx/Ky;

		Fmu(:,i) = dA*(y-mean(y))+mean(y);
		Fsd(:,i) = sqrt(diag(Dss - dA*Dsx'));
		
		if nargout > 2
			Fcov{i} = Dss - dA*Dsx';
		end
	end
end


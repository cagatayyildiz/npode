% compute matrix-valued kernel
%
% uses 'gp.kernel' to determine which kernel to compute:
% - decomposable
% - div-free
% - curl-free
% - div+curl-free
%
function K = matrixkernel(gp, X1, X2, tol, deriv)

	if ~exist('X1','var')
		X1 = gp.X;
	end
	if ~exist('X2','var')
		X2 = gp.X;
	end
	if ~exist('deriv','var')
		deriv = 0;
	end
	if ~exist('tol','var')
		tol = 0;
	end

	
	switch gp.kernel
		case 'gauss'
			K = gausskernel(X1,X2,gp.ell,gp.sf,tol,deriv);
		case 'ma52'
			K = matern52kernel(X1,X2,gp.ell,gp.sf,tol);
		case 'exp'
			K = expkernel(X1,X2,gp.ell,gp.sf,tol,deriv);
		case 'decexp'
			K = decexpkernel(X1,X2,gp.A,gp.ell,gp.sf,tol,deriv);
		case 'decgauss'
			K = decgausskernel(X1,X2,gp.A,gp.ell,gp.sf,tol,deriv);
		case 'div'
			K = divkernel(X1,X2,gp.ell,gp.sf,tol);
		case 'curl'
			K = curlkernel(X1,X2,gp.ell,gp.sf,tol);
		case 'divcurl'
			K = divcurlkernel(X1,X2,gp.ell,gp.sf,tol);
end


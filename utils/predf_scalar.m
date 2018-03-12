%
% ODE function optimised for a single point 'x'
%
% for a general function allowing multiple points, use 'predf'
%
function [Fmu,Fcov] = predf_scalar(gp,Xs, Fs)

	if ~exist('Fs','var')
		Fs = gp.F;
	end

	[Nf,D] = size(gp.F);
	Ns = size(Xs,1);

	switch gp.kernel
		case {'exp','gauss','ma52'}
			Ksx = matrixkernel(gp,Xs,gp.X,0); % kernel bw Xs and all ind. points
%			B = Ksx*gp.invKf;
			B = Ksx/gp.Kf;
			
			% loop version
% 			Fmu = zeros(Ns,D);
% 			for i=1:Ns
% 				Fmu(i,:) = B(i,:)*Fs;
% 			end
            Fmu = B*Fs;
			
			if nargout > 1
				Fcov = gp.sf^2 - sum(B.*Ksx,2); % equal to  diag(B*Ksx')
			end
			
			
%			X = gp.X;
%			Fmu = zeros(Ns,D);
%			Fcov = zeros(Ns,D);
%			for j=1:D
%				lj = gp.ell(j);
%				
%				Kf = gp.sf^2 * exp(-0.5* pdist2(X,X).^2 / lj^2) + 0.001*eye(15);
%				Ksx = gp.sf^2 * exp(-0.5* pdist2(Xs,X).^2 / lj^2); 	
%				Kss = gp.sf^2 * exp(-0.5* pdist2(Xs,Xs).^2 / lj^2); 	
%				B = Ksx/Kf;
%				for i=1:Ns
%					Fmu(i,j) = B*Fs(:,j,i);
%				end
%				Fcov = Kss - B*Ksx';
%			end
			
		otherwise
			D = size(Xs,2);

			Ksx = matrixkernel(gp,Xs,gp.X,0);
			B = Ksx*gp.invKf;
			Fmu = reshape(B*gp.F(:), Ns, D);

			% compute (D,D) covariance
			if nargout > 1
				Kss = gp.sf^2 * gp.A; % always true
				Fcov = Kss - B*Ksx';
			end
	end
end


% 
% 
% function [Fmu,Fcov] = predf1(X,F, gp, x)
% 
% 	switch gp.kernel
% 		case 'dec'
% 			Ksx = kron(gp.A, gp.sf^2 * exp(-0.5*sum((x-X)'.^2) / gp.ell^2) );
% 			B = Ksx * gp.invKy;
% 			Fmu = (B*F(:))';
% 
% 			if nargout > 1
% 				Kss = gp.sf^2 * gp.A;
% 				Fcov = Kss - B*Ksx';
% 			end
% 		case 'div'
% 			[N,D] = size(X);
% %			Ksx = divkernel(x,X,gp.ell,gp.sf);
% 
% 			Dn = sum((x-X)'.^2);
% 			Kg = gp.sf^2 * exp(-0.5 * Dn / gp.ell^2);
% 			% (N,D) pairwise differences
% 			Diff = x-X;
% 			% (D,D,N) outer products of differences
% 			Kd = bsxfun(@times, permute(Diff, [2 3 1]), permute(Diff, [3 2 1])) / gp.ell^2;
% 			% (D,D,N) multiply gaussian kernel
% 			K3 = Kd .* permute(Kg, [3 1 2]);
% 			K = reshape(K3, D, N*D) + kron((D-1-Dn/gp.ell^2).*Kg, eye(D));
% 			Ksx = col2row(K,D);
% 
% 			B = Ksx*gp.invKy;
% 
% 			Fmu = (B*F(:))';
% 
% 			if nargout > 1
% 				Kss = gp.sf^2 * (D-1)*eye(D);
% 				Fcov = Kss - B*Ksx';
% 			end
% 			
% 		case 'curl'
% 			Ksx = curlkernel(x,X,gp.ell,gp.sf);
% 			B = Ksx * gp.invKy;
% 			Fmu = (B*F(:))';
% 
% 			if nargout > 1
% 				Kss = curlkernel(x,x,gp.ell,gp.sf);
% 				Fcov = Kss - B*Ksx';
% 			end
% 	end
% end
% 

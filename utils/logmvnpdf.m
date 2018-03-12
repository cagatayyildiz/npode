function [val,term1,term2] = logmvnpdf(x,m,S)
% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% x is DxN, mu is Dx1, Sigma is DxD


%	logdet = @(A) 2*sum(log(diag(chol(A))));
%	logdet = @(A) sum(log(abs(eig(A))));
	logdet = @(A) sum(log(abs(svd(A))));
	
	D = size(x,1);
	const = -0.5 * D * log(2*pi);
	
	try 
		term1 = -0.5 * (x-m)' / S * (x-m);
		term2 = const - 0.5 * logdet(S);
	catch
		% non SDP, return inf
		val = -inf;
		return;
	end

	val = term1 + term2;
	
	if isinf(val)
		return;
	end
	
	if rcond(S) < 1e-16
		return;
	end
	
	
%	xc = bsxfun(@minus,x,m);
%	term1 = -0.5 * sum((xc' / S) .* xc', 2); 	
end

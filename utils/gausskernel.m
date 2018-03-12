% stationary (scalar) gaussian kernel
%
function K = gausskernel(X1,X2, ells, sf, sn, deriv)
	
	[N1,D] = size(X1);
	N2 = size(X2,1);
	
	if ~exist('ells','var')
		ells = ones(1,D);
	end
	if ~exist('sf','var')
		sf = 1;
	end
	if ~exist('sn','var')
		sn = 0;
	end
	if ~exist('deriv','var')
		deriv = 0;
	end

	
	% matrix-shape
	if N1 > 1 && N2 > 1
		K = sf^2 * exp(-0.5 * distard(X1,X2,ells).^2);
	% row-shape
	elseif N1 == 1 && N2 > 1
		K = sf^2 * exp(-0.5 * distard(X1,X2,ells).^2 ); % optimise
	% single value
	elseif N1 == 1 && N2 == 1
		K = sf^2;
	end
	
	if N1 == N2
		K = K + sn^2*eye(N1);
	end
	
	if deriv == 1
		dK = zeros(N1,N2,D); % tensor
		for j=1:D
			diff = X1(:,j)-X2(:,j)';
			dK(:,:,j) = -ells^-2 * (diff .* K);
		end
		K = dK;
	end
	if deriv == 2
		K = ells^-4 * ( (ells^2 - X1.^2 - X2'.^2 + 2*X1*X2') .* K);
	end
end




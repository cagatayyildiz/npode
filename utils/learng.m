%
% interpolate GP models
%
function gp = learng(t,Y,l)
	
	[N,D] = size(Y);
	gp.D = D;
	gp.N = N;
	gp.t = t;
	gp.Y = Y;
    
    if exist('l','var') && isscalar(l)
        l = l*ones(gp.D,1);
    end
	
    opts = optimoptions('fmincon','display','off'); 
    
	for i=1:D
		y = Y(:,i);
        if nargin == 2 % l learnt
            mll	= @(p) -logmvnpdf(mean(y), y, gausskernel(t,t,p(1),p(2),0) + p(3)^2*eye(N));

            lb = [0 0 0];
            p0 = [1 1 0.1];
            [p,mlli] = fmincon(mll, p0, [], [], [], [], lb,[],[],opts);

            gp.ell(i) = p(1);
            gp.sf(i) = p(2);
            gp.sn(i) = p(3);
            gp.mlli(i) = mlli;
            
        elseif nargin == 3 % l fixed
            mll	= @(p) -logmvnpdf(mean(y), y, gausskernel(t,t,l(i),p(1),0) + p(2)^2*eye(N));

            lb = [0 0];
            p0 = [1 0.1];
            [p,mlli] = fmincon(mll, p0, [], [], [], [], lb,[],[],opts);

            gp.ell(i) = l(i);
            gp.sf(i) = p(1);
            gp.sn(i) = p(2);
            gp.mlli(i) = mlli;
        end
	end
end


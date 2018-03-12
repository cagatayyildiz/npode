function [f,g,pars] = gplvm(Y,q,pars)
%GPLVM 
%
% INPUT
%       Y - original high dimensional data
%       q - number of latent dimensions
%       pars - (sf,sn,l,X)
%
% OUTPUT
%       f - negative log-posterior
%       g - gradient of f

[N,D0] = size(Y);
% assert(length(pars)==N*q+q+2,'parameter vector dimension incorrect');

sf  = exp(pars(1));
sn  = exp(pars(2));
ell = exp(pars(3:q+2))';
X   = reshape(pars(q+3:end),N,q);

Ky = gausskernel(X,X,ell,sf,sn);
Kx = gausskernel(X,X,ell,sf);
f = -sum(diag(logmvnpdf(Y,0,Ky)));

if nargout > 1
    comp_der = @(dk_dth) ( sum(diag(Y'/Ky*dk_dth/Ky*Y))/2 - D0*trace(Ky\dk_dth)/2 );
    g = zeros(length(pars),1);

    dK_dsf = 2*Kx/sf;
    g(1) = comp_der(dK_dsf) * sf;
    
    dK_dsn = 2*sn*(eye(N));
    g(2) = comp_der(dK_dsn) * sn;
    
    for d = 1:q
        dK_dell = X(:,d)-X(:,d)';
        dK_dell = dK_dell.^2 / ell(d)^3 .* Kx;
        g(d+2) = comp_der(dK_dell) * ell(d);
    end
    
    a = Ky\Y;
    invKy = inv(Ky);
    dK = zeros(N,N,q);
    for d=1:q
        dK(:,:,d) = (X(:,d)-X(:,d)') / ell(q)^2  .* Kx;
    end
    dX = zeros(N,q);
    for i=1:D0
        for d=1:q
            dX(:,d) = dX(:,d) + diag( (a(:,i)*a(:,i)' - invKy) * dK(:,:,d) );
        end
    end
    g(q+3:end) = dX(:);
    g = -g; 
end

end
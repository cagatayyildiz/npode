function [ f,g ] = log_gam( x,k,th )
% GAMM computes log-pdf and log-derivative of gamma distribution
%   x  - ndarray
%        points for pdf evalutation
%   k  - postive real
%        shape parameter
%   th - postive real
%        scale parameter
%
%   f  - log pdf 
%   g  - der. pdf wrt. x

f=0; g=0;
if sum(k)>0 && sum(th)>0
    f = log(gampdf(x,k,th));
    g = (k-1)./x - 1./th;
end
end

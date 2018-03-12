function p = par2vec(gp, optpars)
%PAR2VEC forms a parameter vector based on the order in gp.p_order
%
% INPUT
%       gp - np_de object
%       optpars - optimization parameters
%
% OUTPUT
%       p - parameter vector

if ~exist('optpars','var')
    optpars = gp.optpars;
end

p = [];

for pname = gp.p_order
    if contains(optpars,char(pname))
        v = gp.(char(pname));
        p = [p; v(:)];
    end
end
p = p(:);	
    
end

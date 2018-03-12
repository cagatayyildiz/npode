function gp = vec2par(gp, p)
% VEC2PAR sets the parameters of a np_de object based on the parameter
% vector v and the np_de object
%
% INPUT
%       gp - np_de object
%       p - parameter vector
%
% OUTPUT
%       gp - np_de object with parameters read from p
	
[N,D] = size(gp.X);
Nt = length(gp.Y);

i = 0;
% check np_de_model.p_order for the parameter order
if contains(gp.optpars,'log_sf')
    gp.sf = exp(p(i+1));
    i = i + 1;
end
if contains(gp.optpars,'log_sn')
    gp.sn = exp(p(i+1:i+D)');
    i = i + D;
end
if contains(gp.optpars,'log_ell')
    gp.ell = exp(p(i+1:i+D)');
    i = i + D;
end
if contains(gp.optpars,'x0')
    gp.x0 = reshape(p(i+1:i+Nt*D),Nt,D);
    i = i + Nt*D;
end
if contains(gp.optpars,'X')
    gp.X = reshape(p(i+1:i+N*D),N,D);
    i = i + N*D;
end
if contains(gp.optpars,'Fw')
    gp.Fw = reshape(p(i+1:i+N*D),N,D);
    i = i + N*D;
end
end


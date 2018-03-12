function g = wrap_params(gp,varargin)
% WRAP_PARAMS given string-matrix pairs, forms a parameter vector based on
% the order in the np_de object
%
% INPUT
%       gp - np_de object
%       varargin - string-matrix-string-matrix-string-matrix...
%
% OUTPUT
%       g - the vector of parameters in varargin ordered by gp.p_order

g = []; 

for p_ = gp.p_order
    for v=1:length(varargin)/2
        if contains(varargin{v*2-1},p_)
            g = [g; varargin{v*2}(:)];
        end
    end
end 

end

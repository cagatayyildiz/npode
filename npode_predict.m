function X = npode_predict(gp,t,x0)
% NPODE_PREDICT Predicts future using npODE algorithm
%
% t stores the prediction time points. If there is a single trajectory to
% be predicted, t should be a column vector. In case of multiple
% trajectories t is a cell array of column vectors.
%
% Similarly, the output is either a matrix or a cell array of matrices.
% INPUT 
%       gp      - np_de_model object that stores the parameters that are to
%                 be used for prediction
%       t       - prediction time points
%       x0      - initial values, stored in the rows of x0
%
% OUTPUT
%       X       - cell array of predictions

if ~exist('t','var') 
    t = gp.t;
    x0 = gp.x0;
end
sing_ = 0;
if ~iscell(t) && numel(t)==length(t) 
    sing_ = 1;
    t = {t(:)};
end
assert(length(t)==size(x0,1),"the length of t and the number of rows" + ...
    "of x0 must be the same")
X = np_ode_mean_path(gp,t,x0);
if sing_
    X = X{1};
end
end
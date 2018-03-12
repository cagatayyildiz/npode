function [N,Y,X,u,v,X_tr,Y_tr,segments,obs_idx,miss_idx,Y_pca_rec] ...
    = read_data(fname,exp_type,step,q)
%READ_DATA Reads walking sequence data 
%
%   INPUT
%       fname - name of the amc file to be used in the experiment
%       exp_type - 1 for imputation, 2 for forecasting experiments
%       step - downsampling rate
%       q - number of latent dimensions
%
%   OUTPUT
%       N - number of data points in the data set (training + testing)
%       Y - original (high dim) downsampled data (training + testing)
%       X_tr - low dimensional (latent) representation of the data
%              (training + testing)
%       u - principal components
%       v - eigenvalues associated with principal components
%       X_tr - low dimensional (latent) representation of the training data
%       Y_tr - original training data
%       segments - starting indices of the time series given as input
%       obs_idx - indices of the missing part
%       miss_idx - indices of the missing part
%       Y_pca_rec - PCA reconstruction

N = floor(size(amc_to_matrix(fname),1) / step);

if exp_type == 1
    files = {[fname],[fname]};
    [Y,~,~,~] = loadMocapData({[fname]},1,step,N*step);
    C = floor(N/5);
    starts = [1,floor(N/2+C/2)];
    ends = [floor(N/2-C/2), N];
    segments = [1 floor(N/2-C/2)+1];
elseif exp_type == 2
    files = {[fname]};
    [Y,~,~,segments] = loadMocapData({[fname]},1,step,N*step);
    starts = [1];
    ends = [floor(N/2)];
end
D = size(Y,2);

filled_idx = [];
obs_idx = starts(1):ends(1);
for i = 1:length(segments)-1
    filled_idx = [filled_idx ends(i)+1:starts(i+1)-1];
    obs_idx = [obs_idx starts(i+1):ends(i+1)];
end
future_idx =  ends(length(files))+1:N;
miss_idx = [filled_idx future_idx];

% PCA
% Y_true = Y; 
Y_mean = mean(Y);
Y = Y - repmat(Y_mean, N, 1);
[v, u] = pca(Y);
v(find(v<0))=0;
fprintf('Prop. of first %d evalues is %.3f.\n',q,sum(v(1:q)) / sum(v));
X = Y * u(:, 1:q) * diag(1./sqrt(v(1:q)));
Y_pca_rec = X * diag(sqrt(v(1:q))) * u(:, 1:q)';

% training data
Y_tr = Y(obs_idx,:);
X_tr = X(obs_idx,:);
% testing data
% Y_test = Y(miss_idx,:);
% N_test = N - N_tr;

end


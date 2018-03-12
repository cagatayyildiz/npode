function r = logMeanExp(x)
% LOGMEANEXP Numerically stable computation of log(mean(exp(x)))
% 
%   r = LOGMEANEXP(x) computes the log(mean(exp(x), dim)) in a
%   numerically stable way where x is the input vector.
% 
%   Algorithm 
%     Since the computation of the given formula may result to  numerical
%     underflow, we need to make a trick to evaluate the equation as
%     follows:
%     
%     \log\sum\limits_k e^{x_k} 
%      & = \log\left[\left(\sum\limits_k e^{x_k}\right)e^{-y} e^y \right] \\
%      & = \log\left[\left(\sum\limits_k e^{x_k}\right)e^{-y}\right] + y \\
%      & = \log\left[\sum\limits_k e^{x_k - y}\right] + y
% 
%     Choosing $y=\max\limits_k x_k$ in the equation above seems to solve
%     the problem of underflow.
% 
%   Copyright 2011 Ismail Ari.
%   Initial code provided by Baris Kurt and Ali Taylan Cemgil.


    assert(max(size(x)) == numel(x), 'Please provide a 1D array.')

    x = x(:);
    maxX = max(x);    
    r = maxX + log(mean(exp(x-maxX)));
    
end


                              
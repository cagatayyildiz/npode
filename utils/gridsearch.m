function v = gridsearch(vf_gp,fpost)
%GRIDSEARCH Performs grid search over sf, sn and sc and returns
%   the parameter vector giving the highest posterior
%
% INPUT
%       vf_gp - npDE model
%       fpost - function that computes the posterior
%
% OUTPUT
%       v - parameter vector maximizing posterior

sns = linspace(0.2,1,5);
sfs = [1]; 
scs = linspace(0.3,4,5);
f_best = -Inf;
gp_best = vf_gp;
for i=1:length(sns)
    vf_gp.sn = sns(i)*ones(1,vf_gp.D);
    for j=1:length(sfs)
        vf_gp.sf = sfs(j);
        for d = 1:vf_gp.D
            for k=1:length(scs)
                vf_gp = estimate_init_ind_vec(vf_gp);
                vf_gp.Fw(:,d) = scs(k) * vf_gp.Fw(:,d);
                try
                    f_s = fpost(vf_gp);
                    if f_s > f_best
                        f_best = f_s;
                        gp_best = vf_gp;
                    end
                end
            end
        end
    end
end
v = par2vec(gp_best);
end


classdef np_de_model
% NP_DE_MODEL collects all parameters related to the nonparametric
% differential equation solver
%
%  - data (t,Y)
%  - vector field (X,F)
%  - kernel parameter (ell,sf,sn, etc.)
%  - auxiliary variables (optimisation choice, etc.)

	properties
        Nt  % int  - number of time series
		Ny  % vec  - number of data points in each time series - Nt x 1
		x0  % matr - initial points - Nt vec. of size Nt x D
		t   % cell - time points in each time series - Nt vec. of size Ny{i} x 1
        Y   % cell - observations  - Nt matr. of size Ny{i} x D
        
        D   % int  - dimensionality
        Nf  % int  - number of ind. points
		
		% vectors
        W   % int  - grid width
		X   % matr - inducing points - W^D x D
        Fw  % matr - white ind. vecs - W^D x D
		
		% parameters
		log_ell % 1 x D
        log_sf  % 1 
        log_sf0 % initial point 
        log_sn  % 1 x D
        log_sn0 % initial point 
		Kf      % kernel between ind. points
        invKf
        L       % chol(gp.Kf)'
		tol = 0.01;
        
        % hyper-priors
        sf_alpha;  % real
        sf_beta;   % real
        sn_alpha;  % vector
        sn_beta;   % vector
        ell_alpha; % vector
        ell_beta;  % vector
        
		% others 
        ode_fun
		kernel = 'gauss';
		optpars = 'log_sf-log_sn-x0-Fw';
        Xlocs = 'rect';
        p_order = ["log_sf","log_sn","log_ell","x0","X","Fw"];
        clrs = ['r','g','b','m','k','c','y','r','g'];
        non_cent = true;
        x0_true;
        lp = -Inf;
        seed;
    end
    
	% properties to be computed on the fly
	properties (Dependent)
		F
		ell,sf,sn
	end
	methods
		% constructor
		function gp = np_de_model(t,Y,opts)
            
            fnames = fieldnames(opts);
            is_opt_set = @(opt) (~isempty(find(contains(fnames,opt),1)));
            
            gp = gp.set_data(t,Y);
			
            % default values for fields 
            if is_opt_set('W') 
                gp.W = opts.W;
            else
                gp.W = 5;
            end
			gp.Nf = gp.W^gp.D;
            gp.sf_alpha = 0;
            gp.sf_beta = 0;
            gp.sn_alpha = zeros(1,gp.D);
            gp.sn_beta = zeros(1,gp.D);
            gp.ell_alpha = zeros(1,gp.D);
            gp.ell_beta = zeros(1,gp.D);
			gp.log_ell = log( ones(1,gp.D) );
			gp.log_sf = log(1);
            Y_ = [];
            for i=1:gp.Nt, Y_ = [Y_;Y{i}]; end
			gp.log_sn = log((max(Y_) - min(Y_)) / 20);
            
            % inducing point locations
            if is_opt_set('Xlocs')
                gp.Xlocs = opts.Xlocs;
            end
            gp = gen_ind_point_locs(gp);
            
            % setting fields given as argument
            for n = 1:length(fnames)
                try
                    gp.(fnames{n}) = opts.(fnames{n});
                catch
                    fprintf('No np_de_model field for the option %s\n',fnames{n});
                end
            end

            % computing the kernel
			gp.Kf = matrixkernel(gp, gp.X, gp.X, gp.tol);
			gp.invKf = inv(gp.Kf);
			gp.L = chol(gp.Kf)';
            
            % step-0: initial inducing vectors
            gp = gp.estimate_init_ind_vec();
            
            if ~is_opt_set('sf')
                gp.sf = sum(abs(gp.F(:))) / numel(gp.F);
%                 gp.sf = max(abs(gp.F(:)));
            end
            
            gp.log_sf0 = gp.log_sf;
            gp.log_sn0 = gp.log_sn;
        end
        
        function gp = set_data(gp,t,Y)
			gp.t = t;
			gp.Y = Y;
            gp.Nt = length(gp.Y); 
			gp.D = size(gp.Y{1},2);
            gp.Ny = zeros(gp.Nt,1);
            for i = 1:gp.Nt
                gp.Ny(i) = size(gp.Y{i},1);
            end
            if isempty(gp.x0)
                gp.x0 = zeros(gp.Nt,gp.D);
                for i = 1:gp.Nt
                    gp.x0(i,:) = gp.Y{i}(1,:);
                end
            end
		end
        
        
        function gp = gen_ind_point_locs(gp)
            if strcmp(gp.Xlocs,'minbb')  % bounding box
                if gp.D == 2
                    Z = [];
                    for i=1:gp.Nt, Z = [Z;gp.Y{i}]; end
                    Z = Z';
                    bb = minBoundingBox(Z);
                    d1 = bb(:,2)-bb(:,1);
                    l1 = repmat(d1/(gp.W-1),1,gp.W) .* repmat(0:gp.W-1,gp.D,1);
                    l1 = reshape(repmat(l1',1,gp.W)',gp.D,gp.W^2);
                    d2 = bb(:,3)-bb(:,2);
                    l2 = repmat(d2/(gp.W-1),1,gp.W) .* repmat(0:gp.W-1,gp.D,1);
                    l2 = reshape(repmat(l2,1,gp.W),gp.D,gp.W^2);
                    gp.X = l1 + l2 + bb(:,1);
                    gp.X = gp.X';
                else
                    msgID = 'VFMODEL:gen_ind_point_locs';
                    msg = 'min bounding box can only be used in 2D';
                    baseException = MException(msgID,msg);
                    throw(baseException)
                end
            elseif strcmp(gp.Xlocs,'traj')
                gp.Nf = gp.W;
                idx = round(linspace(1,gp.Ny,gp.Nf));
                gp.X = gp.Y{1}(idx,:);
                gp.X = gp.X + randn(size(gp.X))*gp.Nf/500;
            elseif strcmp(gp.Xlocs,'rect') || strcmp(gp.Xlocs,'grid')
                xs = zeros(gp.W,gp.D);
                tmp = [];
                for i=1:gp.Nt, tmp = [tmp;gp.Y{i}]; end
                for i=1:gp.D
                    xs(:,i) = linspace(min(tmp(:,i)),max(tmp(:,i)),gp.W)';
                end
                if gp.D == 2
                    [X1,X2] = meshgrid(xs(:,1), xs(:,2));
                    gp.X = [X1(:) X2(:)];
                elseif gp.D == 3
                    [X1,X2,X3] = meshgrid(xs(:,1), xs(:,2), xs(:,3));
                    gp.X = [X1(:) X2(:) X3(:)];
                end
            end
		end
        
        
        function gp = estimate_init_ind_vec(gp)
            % empirical derivatives at the locations
            F_  = zeros(gp.Nf,gp.D);
            for j = 1:gp.Nt
                dY = diff(gp.Y{j}) ./ diff(gp.t{j});
                Kyy = gausskernel(gp.Y{j}(1:end-1,:),gp.Y{j}(1:end-1,:),gp.ell,gp.sf,mean(gp.sn));
                Y_ = gp.Y{j}(1:end-1,:);
                K = gausskernel(gp.X,Y_,gp.ell,gp.sf,mean(gp.sn));
                F_ = F_ + K / Kyy * dY;
            end
%             F_  = 0.2*rand(gp.Nf,gp.D);
            gp.Fw = gp.L \ F_ / gp.Nt;
		end
		
		% all getters that do something
		function F = get.F(gp)
			switch gp.kernel
				case {'exp','ma32','ma52','gauss'}
					F = gp.L * gp.Fw;  
				case {'div','curl','divcurl','decexp'}
					F = reshape(gp.L * gp.Fw(:),[],gp.D);
			end
        end
		function gp = set.F(gp,F)
			switch gp.kernel
				case {'exp','ma32','ma52','gauss'}
					gp.Fw = gp.L \ F;
				case {'div','curl','divcurl','decexp'}
					gp.Fw = reshape(gp.L\F(:), [], gp.D);
			end
        end

		function gp = update_kernel(gp)
%             if ~gp.non_cent, oldF=gp.F; end
% 			gp.Kf = matrixkernel(gp, gp.X, gp.X, gp.tol);
			gp.Kf = matrixkernel(gp, gp.X, gp.X, 1e-6);
			gp.invKf = inv(gp.Kf);
            try
                gp.L = chol(gp.Kf)';
            catch
                msgID = 'VFMODEL:update_kernel';
                msg = sprintf('Cholesky decomposition cannot be computed. Check gp:\n%s',to_str(gp));
                baseException = MException(msgID,msg);
                throw(baseException)
            end
%             if ~gp.non_cent, gp.F=oldF; end
        end
        
        
        % gp.sf = 1;    gp.('sf') = 1;
		function sf  = get.sf(gp),  sf  = exp(gp.log_sf);  end
		function sn  = get.sn(gp),  sn  = exp(gp.log_sn);  end
		function ell = get.ell(gp), ell = exp(gp.log_ell); end
		
        
		function gp = set.sn(gp,sn)
            if sn < 0, warning('sn < 0'), end
            gp.log_sn = log(sn);
        end 
		function gp = set.sf(gp,sf)
            if sf < 0, warning('sf < 0'), end
            gp.log_sf = log(sf);
            gp = gp.update_kernel();
        end 
		function gp = set.ell(gp,ell) 
            if any(ell < 0), warning('ell < 0'), end
            gp.log_ell = log(ell); 
            gp = gp.update_kernel();
        end 
        
        function str = to_str(gp)
            if gp.D == 2
                str = sprintf('[sf %.2f]\n[ell %.3f %.3f]\n[sn %.3f %.3f]\n',...
                gp.sf,gp.ell(1), gp.ell(2), gp.sn(1), gp.sn(2));
            elseif gp.D == 3
                str = sprintf('[sf %.2f]\n[ell %.3f %.3f %.3f]\n[sn %.3f %.3f %.3f]\n',...
                gp.sf,gp.ell(1), gp.ell(2), gp.ell(3), gp.sn(1), gp.sn(2), gp.sn(3));
            end
        end
        
        function [] = disp(gp)
            fprintf(to_str(gp));
        end
	end
end





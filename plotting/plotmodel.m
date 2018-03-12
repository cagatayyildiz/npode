% plot the ODE system
%
function [] = plotmodel(gp,x0,truef,Nsamples)
		
    if ~exist('x0','var')
		x0 = [];
    end
    if ~exist('Nsamples','var')
		Nsamples = 0;
    end
    if ~exist('truef','var')
		truef = [];
    end
    
    if gp.D == 2
        if gp.Nt == 1
    		subplot(121);
            plotvf(gp,x0,truef);
            subplot(122);
            plotode(gp,1,x0,truef,Nsamples);
        else
            subplot(2,4,[1 2 5 6]);
            plotvf(gp,x0,truef);
            idx = [3,4,7,8];
            for d=1:min(gp.Nt,4)
                subplot(2,4,idx(d));
                plotode(gp,d,x0,truef,Nsamples);
            end
        end
        
    else
        drawnow;
        fs = min(gp.Nt,4);
        w = ceil(sqrt(fs-1))+1;
        for d=1:fs
            subplot(w,w,d);
            plotode(gp,d,x0,truef,Nsamples);
        end
    end
    
end



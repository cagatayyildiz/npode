function [Y, initY, varY, segments] = loadMocapData(files, startFrames, skipFrames, endFrames, expMap)
% implementation in GPDM library: http://www.dgp.toronto.edu/~jmwang/gpdm/

if (~exist('expMap','var'))
    expMap = 0;
end

if (expMap == 0)
    dims =[6 3 3 3 3 3 3 2 3 1 1 2 1 2 2 3 1 1 2 1 2 3 1 2 1 3 1 2 1];
    locations = [1 7 10 13 16 19 22 25 27 30 31 32 34 35 37 39 42 43 44 46 47 49 52 53 55 56 59 60 62];
    mask = [32:33 34 35:36 44:45 46 47:48 55 62];
    %mask = [];
else
    dims =[6 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];
    locations = [1 7 10 13 16 19 22 25 28 31 34 37 40 43 46 49 52 55 58 61 64 67 70 73 76 79 82 85 88];
    mask = [1:6 37:39 40:42 43:45 58:60 61:63 64:66 76:78 88:90];
end



N = size(files,2);
Y = [];
curY = {};
end_segments = [];
for n = 1:N
    if (expMap == 0)
        curY{n} = amc_to_matrix(cell2mat(files(n)));
    else
        curY{n} = amc_to_matrix_exp(cell2mat(files(n)));
    end
    curY{n}(:,mask) = [];
    
  
    Y = [Y; curY{n}];
    if (n == 1)
        end_segments = size(Y,1);
    else
        end_segments = [end_segments end_segments(end) + size(curY{n},1)];
    end
end
initY = Y(1,:);

globY = Y(2:end, 1:3) - Y(1:end-1, 1:3);
globY(end_segments,1:3) = globY(end_segments-1,1:3);
Y = [globY Y(:, 4:end)];

varY = zeros(1,2);
numJoints = (size(Y,2) - 3);
varY(1) = mean(var(Y(:,7:end)));
varY(2) = mean(var(Y(:,4:6)));
varY(3) = mean(var(Y(:,1:3)));


Y = [];
for n = 1:N
    if (endFrames(n) > 0)
        curY{n} = curY{n}(startFrames(n):skipFrames(n):endFrames(n), :);
    else
        curY{n} = curY{n}(startFrames(n):skipFrames(n):end, :);
    end
    
    Y = [Y; curY{n}];
    if (n == 1)
        end_segments = size(curY{n},1);
    else
        end_segments = [end_segments end_segments(end) + size(curY{n},1)];
    end
end

globY = Y(2:end, 1:3) - Y(1:end-1, 1:3);
globY(end_segments,1:3) = globY(end_segments-1,1:3);
Y = [globY Y(:, 4:end)];

% choose one of 2 rescaling schemes

% %next 2 if statements give you rescaling only if thresholds exceeded
% if (varY(1) < varY(2))
%     Y(:,4:6) = Y(:,4:6)*sqrt(varY(1)/varY(2));
% end
% if (varY(1) < varY(3))
%     Y(:,1:3) = Y(:,1:3)*sqrt(varY(1)/varY(3));
% end

% this rescales always 
% Y(:,4:6) = Y(:,4:6)*sqrt(varY(1)/varY(2));
% Y(:,1:3) = Y(:,1:3)*sqrt(varY(1)/varY(3));

segments = end_segments + 1;
segments = [1 segments];
segments(end) = [];


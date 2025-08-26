function cmap = cmapTrans(X,opt)
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
arguments (Repeating)
    X  %(1,:) double {mustBeInRange(X,0,1)} 
end
arguments
    opt.Resolution = 256;
end
n = numel(X);
if n == 1
    cmap = repmat(X{1},opt.Resolution,1);
    return
end
transRange  = floor((opt.Resolution-1)/(n-1))+1;
transA      = linspace(1,0,transRange)';
transB      = linspace(0,1,transRange)';
cmap        = zeros(opt.Resolution,3);
midx        = round(linspace(1,opt.Resolution,n));
for k = 1:n
    if ~isnumeric(X{k})
        switch X{k}
            case {'r', 'red'}
                X{k} = [1 0 0];
            case {'g', 'green'}
                X{k} = [0 1 0];
            case {'b', 'blue'}
                X{k} = [0 0 1];
            case {'k', 'black'}
                X{k} = [0 0 0];
            case {'w', 'white'}
                X{k} = [1 1 1];
            case {'c', 'cyan'}
                X{k} = [0 1 1];
            case {'m', 'magenta'}
                X{k} = [1 0 1];
            case {'y', 'yellow'}
                X{k} = [1 1 0];
        end
    end
    if k < n
        cmap(midx(k):midx(k)+transRange-1,:) = cmap(midx(k):midx(k)+transRange-1,:)+transA*X{k};
        if k > 1
            cmap(midx(k)-transRange+1:midx(k)-1,:) = cmap(midx(k)-transRange+1:midx(k)-1,:)+transB(1:end-1)*X{k};
        end
    else
        cmap(midx(k)-transRange+1:midx(k),:) = cmap(midx(k)-transRange+1:midx(k),:)+transB*X{k};
    end
end
end
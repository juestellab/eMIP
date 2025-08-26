function im = imagescCMAP(imk, colormapval, y, x)
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
if nargin == 4 && all(size(y)>1)
    yTemp = imk;
    xTemp = colormapval;
    imk = y;
    colormapval = x;
    y = yTemp;
    x = xTemp;
end
switch colormapval
    case {'red', 'r'}
        colormapval = [1 0 0];
    case {'green', 'g'}
        colormapval = [0 1 0];
    case {'blue', 'b'}
        colormapval = [0 0 1];
    case {'gray', 'k', 'w'}
        colormapval = [1 1 1];
end

if isnumeric(colormapval) && any(size(colormapval)<=1)
    if size(colormapval,2) == 3 && all(colormapval>=0) && all(colormapval<=1)
        v = (0:255)'./255;
        colormapval = v.*colormapval;
    else       
        rgb = zeros([size(imk),3]);
        rgb(:,:, colormapval) = repmat(imk,1,1,numel(colormapval));
        imk = rescale(rgb);
    end
end

if nargin == 4
    im = imagesc(y,x, imk);
else
    im = imagesc(imk);
end

if ~isnumeric(colormapval) || (size(colormapval,2)==3 && size(colormapval,1)>1)
    colormap(colormapval)
end
end


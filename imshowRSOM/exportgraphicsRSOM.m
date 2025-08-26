function exportgraphicsRSOM(FilePath, FigureNum, Format, opt)
%EXPORTGRAPHICSRSOM 
%   Shortcut to exprot figures with some default settings. For more
%   elaborate export whishes use exportgraphices directly.
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
arguments
    FilePath
    FigureNum                                       
    Format     {mustBeMember(Format, {'.png', '.pdf', '.jpg', '.jpeg', '.tif', '.tiff', '.gif', '.emf', '.eps'})}   = {'.png'}; % example
    opt.ExtensionName                                                                                               = {'_z', '_x', '_y'};
    opt.Resolution                                                                                                  = 300;
    opt.BackgroundColor                                                                                             = 'current';
end

if iscell(Format)
   J = numel(Format);
else
    J = 1;
    Format = {Format};
end

for j = 1:J
    for k = 1:numel(FigureNum)
        set(0, "CurrentFigure", FigureNum(k)); 
        f = gcf;
        exportgraphics(f, [char(FilePath) char(opt.ExtensionName{k}) char(Format{j})], 'Resolution', opt.Resolution, 'BackgroundColor', opt.BackgroundColor);
    end
end
end

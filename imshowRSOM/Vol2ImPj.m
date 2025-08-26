function Im = Vol2ImPj(V, opt)
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
arguments
    V                                                                                               = [];
    opt.axisDir         (1,1) double {mustBeMember(opt.axisDir,[1,2,3])}                            = 1;
    opt.ProjectionType        char {mustBeMember(opt.ProjectionType, {'MIP', 'Energy', 'Sum'})}     = 'MIP';
end
    switch opt.ProjectionType
    case 'MIP'
        Im = max(V,[],opt.axisDir);
    case 'Energy'
        Im = sum(V.^2,opt.axisDir)/size(V,opt.axisDir);
    case 'Sum'
        Im = sum(abs(V),opt.axisDir)/size(V,opt.axisDir);
    end
    Im = squeeze(Im);
end


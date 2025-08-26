function axisDir = pjDir2axisDir(pjDir, opt)
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
arguments
    pjDir                   = [];
    opt.Type       {mustBeMember(opt.Type, {'R.mat', 'alphabetic'})}         = 'R.mat';
end
switch opt.Type
case 'R.mat'
    switch pjDir
    case 'x'
        axisDir = 2;
    case 'y'
        axisDir = 3;
    case 'z'
        axisDir = 1;
    end
case 'alphabetic'
    switch pjDir
    case 'x'
        axisDir = 1;
    case 'y'
        axisDir = 2;
    case 'z'
        axisDir = 3;
    end
end
end


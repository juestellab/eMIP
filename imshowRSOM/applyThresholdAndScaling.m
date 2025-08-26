function Pj = applyThresholdAndScaling(Pj, th, opt)
    if th(1) ~= 0
        Pj          = Pj-th(1);
        th(2)       = th(2) - th(1);
    end
    Pj(Pj<0)        = 0;
    Pj(Pj>th(2))    = th(2);
    Pj              = Pj/(th(2)/255);
    if nargin >= 3
        switch opt.ContrastScaling
        case 'sqrt'
            Pj          = sqrt(Pj);
            Pj          = Pj/(sqrt(255)/255);
        case 'log'
            Pj          = log10(Pj+1);
            Pj          = Pj/(log10(255+1)/255);
        case 'power'
            Pj          = Pj.^opt.ContrastScalingPower;
            Pj          = Pj/(255.^opt.ContrastScalingPower/255);
        end
    end
end

function opt = parseContrast(opt, N)
    if nargin < 2
        N = 3;
    end
    opt.ContrastAlphaHigh       = contrast2cell(opt.ContrastAlphaHigh, N);
    opt.ContrastAlphaLow        = contrast2cell(opt.ContrastAlphaLow, N);
    opt.ContrastFactorLow       = extendND(opt.ContrastFactorLow, N);
    opt.ContrastFactorHigh      = extendND(opt.ContrastFactorHigh, N); 
    opt.ContrastInputCell= cell(N,10);
    for k = 1:N
        opt.ContrastInputCell(k,:) = optContrast2Input(opt,k);
    end
end


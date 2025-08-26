function th = calcThreshold(Pjz, opt)
arguments
    Pjz 
    opt.ContrastAlphaLow                                {mustBeMemberOrInRange(opt.ContrastAlphaLow, {'clipped'},0,1)}                                              = 'clipped';        % qL = ContrastAlphaLow-quantile  of Z-MIP or 0 (clipped);
    opt.ContrastAlphaHigh                               {mustBeInRange(opt.ContrastAlphaHigh, 0,1)}                                                                 = 0.95;             % qH = ContrastAlphaHigh-quantile of Z-MIP;        
    opt.ContrastFactorMethod                            {mustBeMember(opt.ContrastFactorMethod, {'relative', 'absolute'})}                                          = 'relative';       % method how to multiply factors to get final contrast (saturation) threshold
    opt.ContrastFactorLow                                                                                                                                           = 0;                % relative: th_low = ContrastFactorLow*(qH-qL)+qL;   absolute: th_low = ContrastFactorLow*qH;
    opt.ContrastFactorHigh                                                                                                                                          = 1.25;             % relative: th_high= ContrastFactorHigh*(qH-qL)+qL   absolute: th_high= ContrastFactorHigh*qL
    opt.opt = [];
end
    opt = parseInputOpt(opt);
    if isstruct(Pjz)
        Pjz = Pjz.R;
    end
    if numel(size(Pjz)) > 2
        Pjz = Vol2ImPj(Pjz);
    end

    th(2) =  prctile(Pjz, 100*opt.ContrastAlphaHigh, 'all');
    if ~isnumeric(opt.ContrastAlphaLow) && strcmpi(opt.ContrastAlphaLow, 'clipped')
        th(1) = 0;
    else
        th(1) = prctile(Pjz, 100*opt.ContrastAlphaLow, 'all');
    end
    th(1) = min(th(1), th(2));
    th = saturate(th, opt.ContrastFactorLow, opt.ContrastFactorHigh, opt.ContrastFactorMethod);
end

%% INPUT TESTS
function mustBeMemberOrInRange(a, memberCell, lower, upper, varargin)
    if iscell(a)
        charind = cellfun(@ischar, a);
        try
            b = cell2mat(a(~charind));
        catch
            eidType = 'mustBeMemberOrInRange:notMemberOrInRange';
            msgType = 'Input must be numeric, character vector or cell array containing numeric values or character vectors.';
            throwAsCaller(MException(eidType,msgType))
        end
    else
        charind = 1:length(a);
        b  = a;
    end
    if ~isnumeric(a)
        mustBeMember(a(charind), memberCell)
    end
    if isnumeric(b)
        mustBeInRange(b,lower,upper, varargin{:})
    end
end
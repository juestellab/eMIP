function [AF] = equalizeRSOMRecon(LF,HF, opt)
arguments
    LF  = [];
    HF  = [];
    opt.ContrastAlphaLow                                {mustBeMemberOrInRange(opt.ContrastAlphaLow, {'clipped'},0,1)}                                              = 'clipped';        % qL = ContrastAlphaLow-quantile  of Z-MIP or 0 (clipped);
    opt.ContrastAlphaHigh                               {mustBeInRange(opt.ContrastAlphaHigh, 0,1)}                                                                 = 0.95;             % qH = ContrastAlphaHigh-quantile of Z-MIP;        
    opt.ContrastFactorMethod                            {mustBeMember(opt.ContrastFactorMethod, {'relative', 'absolute'})}                                          = 'relative';       % method how to multiply factors to get final contrast (saturation) threshold
    opt.ContrastFactorLow                                                                                                                                           = 0;                % relative: th_low = ContrastFactorLow*(qH-qL)+qL;   absolute: th_low = ContrastFactorLow*qH;
    opt.ContrastFactorHigh                                                                                                                                          = 1.25;             % relative: th_high= ContrastFactorHigh*(qH-qL)+qL   absolute: th_high= ContrastFactorHigh*qL
    opt.weightHF                                                                                                                                                    = 0.5;
end
AF = [];
if isstruct(LF)
    AF = LF;
    LF = LF.R;
end
if isstruct(HF)
    HF = HF.R;
end
opt = parseContrast(opt,2);

thLF = calcThreshold(LF, opt.ContrastInputCell{1,:});
thHF = calcThreshold(HF, opt.ContrastInputCell{1,:});

R = ((1-opt.weightHF)*applyThresholdAndScaling(LF, thLF) + opt.weightHF*applyThresholdAndScaling(HF, thHF));
if isstruct(AF)
    AF.R = R;
else
    AF = R;
end

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
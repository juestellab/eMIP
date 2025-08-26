function opt = parseReconParams(R, opt)
    if isstruct(R) 
        if isempty(opt.ReconParams)
            if isfield(R, 'reconParams')
                opt.ReconParams = R.reconParams;
            elseif isfield(R, 'ReconParams')
                opt.ReconParams = R.ReconParams;
            elseif isfield(R, 'rP')
                opt.ReconParams = R.rP;
            end
        end
    end
end


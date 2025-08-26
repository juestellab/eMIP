function opt = parseInputOpt(opt, overwrite)
arguments
    opt                                                 = []
    overwrite                                           = true;
end
if ~isfield(opt, 'opt')
    return
end
if ~isstruct(opt.opt)
    opt = rmfield(opt, 'opt');
    return
end
opt2 = opt.opt;
fNames = fieldnames(opt2);
for k = 1:numel(fNames)
    if overwrite || ~isfield(opt, fNames{k})
        opt.(fNames{k}) = opt2.(fNames{k});
    end
end
opt = rmfield(opt, 'opt');
end


function contrast = contrast2cell(contrast,N)
    if nargin < 2
        N = 3;
    end
    if isnumeric(contrast)
        contrast = num2cell(contrast);
    end

    if ~iscell(contrast)
        contrast = {contrast};
    end

    for k = 2:N
        if numel(contrast) < k
            contrast{k} = contrast{k-1};
        end
    end
end


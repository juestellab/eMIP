function [Rshift, surface_pos] = surfaceCorrectionRSOM(R_in, surface, opt)
    % Shifts the RSOM recon according to the surface map and the additional
    % parameter
    % Input:    - RSOM recon
    %           - surface map (2D values)
    %           - parameter:
    %               * 'surface_pos': {'top', 'bottom', interger
    %                   value} (defines where to locate the surface)
    %               * 'image_size': {'max', 'keep_size', integer value}
    %                   (defines how to crop the image)
    % Output:   - shifter Rsom recon (Rshift)
    %           - position of the surface in Rshift 
    arguments
        R_in            = [];
        surface         = [];
        opt.SurfacePos  {mustBeMemberOrInRange(opt.SurfacePos, {'top', 'bottom'}, 0, 1e8)}              = 100;           % position in pixel; top -> locates surface at heighst detected surface point; bottom -> locates surface at lowest detected surface point;
        opt.ImageSize   {mustBeMemberOrInRange(opt.ImageSize, {'keepSize', 'max'}, 1, 1e8)}             = 'keepSize';
    end
    %% parameter
    if isa(surface, 'sfit')
        x = 1:size(R_in,2);
        y = 1:size(R_in,3);
        [X,Y] = meshgrid(x,y);
        X     = X';
        Y     = Y';
        surface = max(1,round(surface(X,Y)));
    end
    
    %% Rshift
    max_surface     = max(surface(:));
    surface_idx     = max_surface - surface;

    Rshift = zeros(size(R_in,1) + max(surface_idx(:)), size(R_in,2),size(R_in,3), 'single');

    for ix=1:size(Rshift,2)
        for iy=1:size(Rshift,3)
            Rshift(surface_idx(ix,iy) + (1:size(R_in,1)),ix,iy) = R_in(:,ix,iy);
        end
    end

    %% Positioning of the surface
    switch opt.SurfacePos
        case 'bottom'
            surface_pos = max_surface+1;
        case 'top'
            surface_pos = min(surface(:))+1;
            max_surface_idx = max(surface_idx(:));
            Rshift = Rshift(1+max_surface_idx:end,:,:);
        otherwise
            surface_pos = opt.SurfacePos;
            idx_correction = max(surface(:)) - surface_pos;
            if idx_correction > 0 
                Rshift = Rshift(idx_correction:end , :, :);
            elseif idx_correction < 0
                Rshift = cat(1,zeros(max(1,-idx_correction), size(Rshift,2), size(Rshift,3)), Rshift);
            end
    end

    %% Cropping the image
    switch opt.ImageSize
        case 'max'
            return
        case 'keepSize'
            if size(R_in,1) > size(Rshift,1)
                 Rshift = cat(1, Rshift, zeros(size(R_in,1)-size(Rshift,1), size(R_in,2), size(R_in,3)));
            else
                Rshift = Rshift(1:size(R_in,1),:,:);
            end
        otherwise
            if opt.ImageSize > size(Rshift,1)
                Rshift = cat(1, Rshift, zeros(opt.ImageSize-size(Rshift,1), size(R_in,2), size(R_in,3)));
            else
                 Rshift = Rshift(1:opt.ImageSize,:,:);
            end               
    end
end

%%
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


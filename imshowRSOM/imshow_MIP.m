function [im1, im_location] = imshow_MIP(R, axis_dir, varargin)
% Some default image adjustment for RSOM recon MIPs
% Input:    - R:     RSOM recon
%           - axis:  direction of MIP
% Varargin: - reconParams: recon parameters (ds, dz)
%           - surface_h_in_pixel/surface_h_in_um: position of surface
%           the axis -> the image location will be adapted such that
%           the surface is at the 0 line
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
    ax = gca;            
    reconParams = [];
    surface_h = 0;
    surface_h_in_pixel = [];
    surface_h_in_mm = [];
    surface_h_in_um = [];
    unit = 'mm';
    daspect_value = [];

    for ind_k = 1:2:numel(varargin)
    switch lower(varargin{ind_k})
        case {'reconparmas', 'reconparams'}
            reconParams = varargin{ind_k+1};
        case 'surface_h_in_pixel'
            surface_h_in_pixel = varargin{ind_k+1};
        case 'surface_h_in_um'
            surface_h_in_um = varargin{ind_k+1};
        case 'surface_h_in_mm'
            surface_h_in_mm =  varargin{ind_k+1};
        case 'unit' 
            unit = varargin{ind_k+1};
        case 'daspect'
            daspect_value = varargin{ind_k+1};
    end
    end

    if ~isempty(reconParams)
        try
            ds = reconParams.GRID_DS;
            dz = reconParams.GRID_DZ;
        catch
            ds = reconParams.ds;
            dz = reconParams.dz;
        end
    else
        ds = 12;
        dz = 3;
    end

    switch unit
    case 'mm'
        ds = ds/1000;
        dz = dz/1000;
        surface_h_in_um = surface_h_in_um/1000;
    case 'um'
        surface_h_in_mm = surface_h_in_mm*1000;
    case 'pixel'
        surface_h_in_um = surface_h_in_um/dz;
        surface_h_in_mm = surface_h_in_mm*1000/dz;
        ds = 1;
        dz = 1;
    end

    surface_h_in_pixel = surface_h_in_pixel*dz;
    
    if ~isempty(surface_h_in_mm)
        surface_h = surface_h_in_mm;
    elseif ~isempty(surface_h_in_um)
        surface_h = surface_h_in_um;
    elseif ~isempty(surface_h_in_pixel)
        surface_h = surface_h_in_pixel;
    end

    if numel(size(R)) == 3
        switch axis_dir
        case 1
            x = [0 ds*size(R,2)];
            y = [0 ds*size(R,3)];
        case 2
            x = [0 dz*size(R,1)]-surface_h;
            y = [0 ds*size(R,3)];
        case 3
            x = [0 dz*size(R,1)]-surface_h;
            y = [0 ds*size(R,2)];
        end
        R = squeeze(max(R,[],axis_dir));
    end

    im_location = [y;x];
    im1 = imagesc(y,x, R);

    if isempty(daspect_value)
        if axis_dir == 1
            daspect([1 1 1])
        else
            daspect([2 1 1])
        end
    elseif isstring(daspect_value) || ischar(daspect_value)
        axis(ax, daspect_value)
    elseif iscell(daspect_value)
        axis(ax, daspect_value{:})
    elseif isnumeric(daspect_value) 
        daspect(daspect_value)
    end
end


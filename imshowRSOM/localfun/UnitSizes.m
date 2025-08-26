function opt = UnitSizes(opt)
    arguments
        opt.ReconParams                         struct                                                                                                                  = [];               % recon parameters (automatically extracted from LF if it is imported structure)
        opt.SurfaceHeight                                   {mustBeNumeric}                                                                                             = [];               % position of the surface -> the image location will be adapted such that the surface is at the 0 line
        opt.SurfaceHeightInPixel                            {mustBeNumeric}                                                                                             = [];               % used for surface height input in pixel
        opt.SurfaceHeightInmm                               {mustBeNumeric}                                                                                             = 0.4;              % used for surface height input in pixel mm
        opt.SurfaceHeightInum                               {mustBeNumeric}                                                                                             = [];               % used for surface height input in pixel um
        opt.Crop                                            {mustBeGreaterThanOrEqual(opt.Crop,0)}                                                                      = [0 0 0; inf inf inf]; % Crop (2,1): [z_start;z_end], (2,2): [z_start,y_start;z_end,y_end] or (2,3):[z_start,y_start,x_start;z_end,y_end,x_end] matrix -> crops 3D-recons accordingly
        opt.CropInPixel                                     {mustBeGreaterThanOrEqual(opt.CropInPixel,0)}                                                               = zeros(0,3);       % Crop given in Pixel unit
        opt.CropInmm                                        {mustBeGreaterThanOrEqual(opt.CropInmm,0)}                                                                  = zeros(0,3);       % Crop given in mm
        opt.CropInum                                        {mustBeGreaterThanOrEqual(opt.CropInum,0)}                                                                  = zeros(0,3);       % Crop given in um
        opt.Unit                                char        {mustBeMember(opt.Unit, {'mm','um', 'pixel'})}                                                              = 'pixel';             % axis Unit
        opt.opt                                                                                                                                                         = [];
    end
    opt = parseInputOpt(opt);
    if ~isempty(opt.ReconParams)
        try
            ds = opt.ReconParams.GRID_DS;
            dz = opt.ReconParams.GRID_DZ;
        catch
            ds = opt.ReconParams.ds;
            dz = opt.ReconParams.dz;
        end
    else
        ds = 12;
        dz = 3;
    end

    switch opt.Unit
    case 'mm'
        ds = ds/1000;
        dz = dz/1000;
        opt.SurfaceHeightInum = opt.SurfaceHeightInum/1000;
        opt.CropInum          = opt.CropInum/1000;
    case 'um'
        opt.SurfaceHeightInmm = opt.SurfaceHeightInmm*1000;
        opt.CropInmm          = opt.CropInmm*1000;
    case 'pixel'
        opt.SurfaceHeightInum = opt.SurfaceHeightInum/dz;
        opt.SurfaceHeightInmm = opt.SurfaceHeightInmm*1000/dz;
        opt.CropInum          = opt.CropInum./[dz ds ds];
        opt.CropInmm          = opt.CropInmm*1000./[dz ds ds];
        ds = 1;
        dz = 1;
    end

    opt.SurfaceHeightInPixel = opt.SurfaceHeightInPixel*dz;
    opt.CropInPixel          = opt.CropInPixel.*[dz ds ds];
    
    if ~isempty(opt.SurfaceHeight)
        ; %#ok<NOSEMI> %do nothing
    elseif ~isempty(opt.SurfaceHeightInPixel)
        opt.SurfaceHeight = opt.SurfaceHeightInPixel;
    elseif ~isempty(opt.SurfaceHeightInum)
        opt.SurfaceHeight = opt.SurfaceHeightInum;
    elseif ~isempty(opt.SurfaceHeightInmm)
        opt.SurfaceHeight = opt.SurfaceHeightInmm;
    end

    if ~isempty(opt.CropInPixel)
        opt.Crop = opt.CropInPixel;
    elseif ~isempty(opt.CropInum)
        opt.Crop = opt.CropInum;
    elseif ~isempty(opt.CropInmm)
        opt.Crop = opt.CropInmm;
    end

    opt.ds = ds;
    opt.dz = dz;
    opt.SurfaceHeightInPixel = round(opt.SurfaceHeight/dz);
    opt.CropInPixel          = round(opt.Crop./[dz ds ds]);
end

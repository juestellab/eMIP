function [R, opt] = cropRecon(R, opt)
arguments
        R = [];
        opt.ReconParams                                     struct                                                                                                                  = [];               % recon parameters (automatically extracted from LF if it is imported structure)
        opt.Crop                                            {mustBeGreaterThanOrEqual(opt.Crop,0)}                                                                      = [0 0 0; inf inf inf]; % Crop (2,1): [z_start;z_end], (2,2): [z_start,y_start;z_end,y_end] or (2,3):[z_start,y_start,x_start;z_end,y_end,x_end] matrix -> crops 3D-recons accordingly
        opt.CropInPixel                                     {mustBeGreaterThanOrEqual(opt.CropInPixel,0)}                                                               = zeros(0,3);       % Crop given in Pixel unit
        opt.CropInmm                                        {mustBeGreaterThanOrEqual(opt.CropInmm,0)}                                                                  = zeros(0,3);       % Crop given in mm
        opt.CropInum                                        {mustBeGreaterThanOrEqual(opt.CropInum,0)}                                                                  = zeros(0,3);       % Crop given in um                opt.Unit                                char        {mustBeMember(opt.Unit, {'mm','um', 'pixel'})}                                                              = 'pixel';             % axis Unit
        opt.Unit                                char        {mustBeMember(opt.Unit, {'mm','um', 'pixel'})}                                                              = 'mm';             % axis Unit
        opt.opt                                                                                                                                                         = [];
        opt.skipUnitSizes                                                                                                                                               = false;                
end


opt = parseInputOpt(opt);
opt = parseReconParams(R, opt);

if isstruct(R) 
    R  = R.R;
end

opt = extendCropInfo(opt);

sz_max = size(R);

if opt.skipUnitSizes
    opt = UnitSizes('opt', opt);     
end

opt.CropInPixel(1,:) = min(max(opt.CropInPixel(1,:),1), sz_max);
opt.CropInPixel(2,:) = min(max(opt.CropInPixel(2,:),opt.CropInPixel(1,:)), sz_max);
pos = opt.CropInPixel;

opt.Crop = [(opt.CropInPixel(1,:)-1).*[opt.dz opt.ds opt.ds]; opt.CropInPixel(2,:).*[opt.dz opt.ds opt.ds]];


if prod(pos(1,:))>1.1 || sum(sz_max-pos(2,:))>0.1
    R = R(pos(1,1):pos(2,1), pos(1,2):pos(2,2), pos(1,3):pos(2,3));
    opt.wasCroped = true;
else
    opt.wasCroped = false;
end

end

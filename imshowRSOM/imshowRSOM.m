function [img, Pj, opt, R_C] = imshowRSOM(LF, HF, opt)
% Visualization for RSOM reconstruction.
% Performs maximum intensity projections along the designated axes ('x',
% 'y' or 'z'). 
% Default visualization:
% Maximum intensity projections (MIPs) along the x-, y- and z-axis for low 
% frequency (LF) and high frequency (HF) reconstruction. The contrast is 
% calculated for LF and HF separately by clipping negative values to 0 and 
% capping values at 1.25 times the 95th percentile value of the respective
% z-MIPs. The corresponding MIPs are then fused into three rgb-images. LF is
% assigned linearly scaled to the red channel and HF to the green channel. 

% Example: 
% img = imshowRSOM(LF, HF, 'Unit', 'um', 'ContrastAlphaHigh', [0.9, 0.92], 
% 'ContrastAlphaLow', {'clipped', 0.025}, 'SurfaceHeightInum', 600)
%
% Output: 
% -img: structure with image objects, e.g.; img.imZ.CData shows the
% 3D-RGB data of the shown z-MIP.
% - Pj: projection data
% - opt: options used to create the plots

% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de

    arguments
        %% Recons
        LF                                                                                                                                                              = [];               % Low freq recon (default red channel) or single recon; data format: recon or recon structure
        HF                                                                                                                                                              = [];               % High freq recon (default green channel), if empty LF -> R; data format: recon or recon structure
        opt.AF                                                                                                                                                          = [];        
        opt.SignalSign                  (1,1)               {mustBeMember(opt.SignalSign,[-1,1])}                                                                       = 1;                % sometimes signal signs are reveresed, with SignalSign = -1 one can undo it.                                               
        opt.ReconParams                         struct                                                                                                                  = [];               % recon parameters (automatically extracted from LF if it is imported structure)
        %% Visualization - Axis
        opt.SurfaceHeight                                   {mustBeNumeric}                                                                                             = [];               % position of the surface -> the image location will be adapted such that the surface is at the 0 line
        opt.SurfaceHeightInPixel                            {mustBeNumeric}                                                                                             = [];               % used for surface height input in pixel
        opt.SurfaceHeightInmm                               {mustBeNumeric}                                                                                             = 0.4;              % used for surface height input in pixel mm
        opt.SurfaceHeightInum                               {mustBeNumeric}                                                                                             = [];               % used for surface height input in pixel um
        opt.Crop                                            {mustBeGreaterThanOrEqual(opt.Crop,0)}                                                                      = [0 0 0; inf inf inf]; % Crop (2,1): [z_start;z_end], (2,2): [z_start,y_start;z_end,y_end] or (2,3):[z_start,y_start,x_start;z_end,y_end,x_end] matrix -> crops 3D-recons accordingly
        opt.CropInPixel                                     {mustBeGreaterThanOrEqual(opt.CropInPixel,0)}                                                               = zeros(0,3);       % Crop given in Pixel unit
        opt.CropInmm                                        {mustBeGreaterThanOrEqual(opt.CropInmm,0)}                                                                  = zeros(0,3);       % Crop given in mm
        opt.CropInum                                        {mustBeGreaterThanOrEqual(opt.CropInum,0)}                                                                  = zeros(0,3);       % Crop given in um
        opt.Unit                                char        {mustBeMember(opt.Unit, {'mm','um', 'pixel'})}                                                              = 'mm';             % axis Unit
        opt.FastScanningAxis                    char        {mustBeMember(opt.FastScanningAxis, {'x', 'y', 'larger'})}                                                  = 'larger';         % defines the fast scanning axis (either x or y), the axis-label will show the fast scanning always as the y-axis
        %% Visaualization - Aspect
        opt.AspectXYZ                                                                                                                                                   = [];               % Aspect ratio [X Y Z]: for Z-MIP - > daspect([X Y 1]), Y-MIP -> daspect([X Z 1], Z-MIP -> daspect([Y Z 1]) 
        opt.DaspectValue                                                                                                                                                = 'default';        % aspect ratio of images 
        opt.DaspectValueX                                                                                                                                               = [];               % overwrites aspect_value for MIP in x-dir 
        opt.DaspectValueY                                                                                                                                               = [];               % overwrites aspect_value for MIP in y-dir 
        opt.DaspectValueZ                                                                                                                                               = [];               % overwrites aspect_value for MIP in z-dir 
        %% Visualization - Contrast
        opt.CMAP                                                                                                                                                        = 'hot';            % colormap for single channel visualization, if numeric -> rgb: red = 1, green = 2, blue = 3, gray = [1 2 3]; 
        opt.ThresholdLF                                                                                                                                                 = [];               % absolute threshold values of LF (All if HF is not defined) for contrast saturation; if defined, contrast saturation and FBE will be bypassed
        opt.ThresholdHF                                                                                                                                                 = [];               % absolute threshold values of HF for contrast saturation; (only active if ThresholdLF is defined)
        opt.ContrastAlphaLow                                {mustBeMemberOrInRange(opt.ContrastAlphaLow, {'clipped'},0,1)}                                              = 'clipped';        % qL = ContrastAlphaLow-quantile  of Z-MIP or 0 (clipped);     scalar, string or 2d cell array (for LF, HF individually)
        opt.ContrastAlphaHigh                               {mustBeInRange(opt.ContrastAlphaHigh, 0,1)}                                                                 = 0.95;             % qH = ContrastAlphaHigh-quantile of Z-MIP;                    scalar         or 2d cell array (for LF, HF individually)
        opt.ContrastFactorMethod                            {mustBeMember(opt.ContrastFactorMethod, {'relative', 'absolute'})}                                          = 'relative';       % method how to multiply factors to get final contrast(saturation) threshold
        opt.ContrastFactorLow                                                                                                                                           = 0;                % relative: th_low = ContrastFactorLow*(qH-qL)+qL;   absolute: th_low = ContrastFactorLow*qH;
        opt.ContrastFactorHigh                                                                                                                                          = 1.25;             % relative: th_high= ContrastFactorHigh*(qH-qL)+qL   absolute: th_high= ContrastFactorHigh*qL
        opt.ContrastScaling                     char        {mustBeMember(opt.ContrastScaling, {'linear', 'sqrt', 'log', 'power'})}                                     = 'linear';         % scaling between th_low and th_high value; 
        opt.ContrastScalingPower        (1,1)               {mustBeNumeric}                                                                                             = 2;                % power factor if scaling is power
        opt.MergeChannels                       logical                                                                                                                 = false;            % Merge LF and HF into one channel (e.g., gray-scale, hot, ...)
        opt.ColorChannels               (1,3)               {mustBeInRange(opt.ColorChannels,0,3), mustBeInteger(opt.ColorChannels)}                                    = [1 2 0];          % RGB coloring channels: 1 = LF, 2 = HF; 
        opt.FBE                         (1,1)   logical                                                                                                                 = false;            % FBE: Frequency Band Equalization, if true: HF = alpha*HF s.t. min_alpha || LF-alpha*HF || 
        opt.FBEoptions                          char        {mustBeMember(opt.FBEoptions, {'Pjx', 'Pjy', 'Pjz', 'Pji', 'volume','energy'})}                             = 'Pjz';            % Pjx_y_z: FBE based on Pj along x,y,z axis; ; Pji: FBE for each Pj individually; volume: FBE based on 3D volume; energy; FBE based on volume squared  
        opt.FBEcontrastMethod                   char        {mustBeMember(opt.FBEcontrastMethod, {'joint', 'independent'})}                                             = 'joint'
        opt.Alphaval                                                                                                                                                    = [];               % predifined Alphaval (only active for FBE and if ThresholdLF is defined)
        %% Surface Correction
        opt.Surface                                                                                                                                                     = [];
        opt.SurfaceVarargin                     cell                                                                                                                    = {};
        %% Figure Properties
        opt.BackgroundColor                                                                                                                                             = 'k';              % Bacground of figure; black = 'k', white = 'w'
        opt.LabelColor                                                                                                                                                  = 'w';              % Color of labels;  black = 'k', white = 'w' 
        opt.FigureNum                                       {mustBeInteger(opt.FigureNum)}                                                                              = [];               % figure number for imZ, imY, imX
        opt.WindowState                                     {mustBeMember(opt.WindowState,{'normal', 'fullscreen', 'minimized', 'maximized'})}                          = 'normal';                                                                                                   
        opt.FigurePosition                      double                                                                                                                  = [];
        opt.FigureUnits                                     {mustBeMember(opt.FigureUnits, {'centimeters', 'characters', 'inches', 'normalized', 'pixels', 'points'})}  = 'pixels';
        opt.Visible                             char        {mustBeMember(opt.Visible, {'on', 'off'})}                                                                  = 'on';             % if off, figures become invisible        
        opt.ProjectionDir                       char        {mustBeMember(opt.ProjectionDir, {'z','x','y','zx','zy','xy','zxy'})}                                       = 'zxy';            % images/MIPS to be visualized
        opt.ProjectionType                      char        {mustBeMember(opt.ProjectionType, {'MIP', 'Energy', 'Sum'})}                                                = 'MIP';            % type of image projection 
        opt.AxisLabelAndTicks                   logical                                                                                                                 = true;
    end
    
    %% Main
    % Initialization
    [LF, HF, AF, opt, img, Pj] = initializeFun(LF, HF, opt);

    % Predfiend contrast thresholds
    [opt, Pj] = predefinedContrast(opt, Pj);

    % Extract Resolution & set units
    opt = UnitSizes(opt = opt);

    % Surface correction
    [LF, HF, AF,  opt] = surfaceCorrection_internal(LF, HF, AF, opt);
    
    % Pj along z direction
    Pj = projection(LF, HF, AF, opt, Pj, 'z'); 

    % Contrast (threshold) calculation (and FBE)
    Pj = contrastCalc(LF, HF, AF, opt, Pj);
    
    % Crop
    [LF, HF, AF, opt, Pj, R_C] = cropRecon_internal(LF, HF, AF, opt, Pj);
   

    % Imaging
    [img, Pj] = imaging_internal(LF, HF, AF, opt, Pj, img);
    
    % Reset default figure values
    set(groot, 'DefaultFigureVisible', opt.defaultFigVis);
end

%% SUPPORT FUNCTIONS
%% Initialization
function [LF, HF, AF, opt, img, Pj] = initializeFun(LF, HF, opt)
    opt.defaultFigVis = get(groot, 'DefaultFigureVisible');
    set(groot, 'DefaultFigureVisible', 'off');
    img                         = [];
    Pj                          = struct('LF',[],'HF',[], 'AF', []);
    opt                         = parseContrast(opt);
    opt.singleFB                = isempty(HF);
    opt.FBE                     = opt.FBE && ~opt.singleFB;
    opt.calcContrast            = true;
    opt.figureIter              = 0;
    AF                          = opt.AF;
    
    opt = extendCropInfo(opt);

    % Structure to 3D Volume & optional sign switch
    [LF, HF, AF, opt] = EvalInitStruct(LF, HF, AF, opt);

    % Fast scanning axis
    if strcmp(opt.FastScanningAxis, 'larger') && size(LF,2) >= size(LF,3)
        opt.FastScanningAxis = 'x';
    elseif strcmp(opt.FastScanningAxis, 'larger')
        opt.FastScanningAxis = 'y';
    end
end

function [LF, HF, AF, opt] = EvalInitStruct(LF, HF, AF, opt)
    opt = parseReconParams(LF, opt);
    if isstruct(LF) 
        LF  = opt.SignalSign * LF.R;
    else
        LF  = opt.SignalSign * LF;
    end
    if isstruct(HF)
        HF  = opt.SignalSign * HF.R;
    else
        HF  = opt.SignalSign * HF;
    end
    if isstruct(AF)
        AF  = opt.SignalSign * AF.R;
    else
        AF  = opt.SignalSign * AF;
    end
end

%% Predfiend contrast thresholds
function [opt, Pj] = predefinedContrast(opt, Pj)
    if ~isempty(opt.ThresholdLF)
        Pj.LF_th = opt.ThresholdLF;
        if ~isempty(opt.ThresholdHF)
            Pj.HF_th = opt.ThresholdHF;
        else
            Pj.HF_th = opt.ThresholdLF;
        end
        if ~isempty(opt.Alphaval)
            Pj.Alphaval = opt.Alphaval;
        else
            Pj.Alphaval = 1;
        end
        opt.calcContrast = false;
    end
end 

%% Surface correction
function [LF, HF, AF, opt] = surfaceCorrection_internal(LF, HF, AF, opt)
    if ~isempty(opt.Surface)
        [LF, opt.SurfaceHeightInPixel] = surfaceCorrectionRSOM(LF, opt.Surface, 'SurfacePos', opt.SurfaceHeightInPixel, opt.SurfaceVarargin{:});
        if ~opt.singleFB
            HF = surfaceCorrectionRSOM(HF, opt.Surface, 'SurfacePos', opt.SurfaceHeightInPixel, opt.SurfaceVarargin{:});
            if ~isempty(AF)
                AF = surfaceCorrectionRSOM(AF, opt.Surface, 'SurfacePos', opt.SurfaceHeightInPixel, opt.SurfaceVarargin{:});
            end
        end
    end
end

%% Crop
function [LF, HF, AF, opt, Pj, R_C] = cropRecon_internal(LF, HF, AF, opt, Pj)
    [LF, opt] = cropRecon(LF, opt = opt);
    if ~opt.singleFB
        [HF, opt] = cropRecon(HF, 'opt', opt, 'skipUnitSizes', true);
        if ~isempty(AF)
            [AF, opt] = cropRecon(AF, opt = opt);
        end
    end
    if opt.wasCroped
        Pj.LF = [];
        Pj.HF = [];
        Pj.AF = [];
    end
    R_C.LF = LF;
    R_C.HF = HF;
    R_C.AF = AF;
end

%% Projection
function Pj = projection(LF, HF, AF, opt, Pj, pjDir)
    axisDir = pjDir2axisDir(pjDir);
    if ~isfield(Pj.LF, pjDir)
        Pj.LF.(pjDir) = Vol2ImPj(LF, "axisDir", axisDir,"ProjectionType", opt.ProjectionType);
    end
    if ~opt.singleFB
        if ~isfield(Pj.HF, pjDir)
            Pj.HF.(pjDir) = Vol2ImPj(HF, "axisDir", axisDir,"ProjectionType", opt.ProjectionType);
        end
        if ~isfield(Pj.AF, pjDir) & ~isempty(AF)
            Pj.AF.(pjDir) = Vol2ImPj(AF, "axisDir", axisDir,"ProjectionType", opt.ProjectionType);
        end
    end
end

%% Contrast & Saturation



function Pj = contrastCalc(LF, HF, AF, opt, Pj)
    if ~opt.calcContrast
        return
    end
    
    switch opt.FBE
    % Seperate Frequancy band Contrast & Saturation
    case false
        Pj.Alphaval = 1;
        Pj.LF_th = calcThreshold(Pj.LF.z, opt.ContrastInputCell{1,:});
        if ~opt.singleFB 
            Pj.HF_th = calcThreshold(Pj.HF.z,opt.ContrastInputCell{2,:});
            if ~isempty(AF)
                Pj.AF_th = calcThreshold(Pj.AF.z, opt.ContrastInputCell{3,:});
            end
        end

    case true
        Pj = FBEcalculation(opt, Pj, LF, HF);
    end
end

function Pj = FBEcalculation(opt, Pj, LF, HF)
% Frequency band equalization
    if ~strcmpi(opt.FBEoptions, 'pji')
    % Frequency band equalization; common thresholds for LF and HF are selected         
        [Alphaval, Pj] = FBE_alpha(LF, HF, 'FBEoptions', opt.FBEoptions, 'Pj', Pj);
        Pj.Alphaval = Alphaval;
        if strcmp(opt.FBEcontrastMethod, 'joint')
            th          = calcThreshold([Pj.LF.z, Pj.Alphaval*Pj.HF.z], 'ContrastAlphaLow', opt.ContrastAlphaLow{1}, 'ContrastAlphaHigh', opt.ContrastAlphaHigh{1}, 'ContrastFactorLow', 1, 'ContrastFactorHigh', 1, 'ContrastFactorMethod', 'absolute');
            Pj.LF_th    = saturate(th, opt.ContrastFactorLow(1), opt.ContrastFactorHigh(1), opt.ContrastFactorMethod);
            Pj.HF_th    = saturate(th, opt.ContrastFactorLow(2), opt.ContrastFactorHigh(2), opt.ContrastFactorMethod);
        else
            Pj.LF_th    = calcThreshold(Pj.LF.z, opt.ContrastInputCell{1,:});
            Pj.HF_th    = calcThreshold(Pj.Alphaval*Pj.HF.z, opt.ContrastInputCell{2,:});
        end
        return
    end

    for k = opt.ProjectionDir
        pjDir                   = lower(k);
        [Alphaval, Pj]          = FBE_alpha(LF, HF, 'FBEoptions', ['Pj' pjDir], 'Pj', Pj);
        Pj.Alphaval.(pjDir)     = Alphaval;
        if strcmp(opt.FBEcontrastMethod, 'joint')
            th                  = calcThreshold([getFieldVal(Pj.LF,pjDir), getFieldVal(Pj.Alphaval,pjDir)*getFieldVal(Pj.HF,pjDir)],  'ContrastAlphaLow', opt.ContrastAlphaLow{1}, 'ContrastAlphaHigh', opt.ContrastAlphaHigh{1}, 'ContrastFactorLow', 1, 'ContrastFactorHigh', 1, 'ContrastFactorMethod', 'absolute');
            Pj.LF_th.(pjDir)    = saturate(th, opt.ContrastFactorLow(1), opt.ContrastFactorHigh(1), opt.ContrastFactorMethod);
            Pj.HF_th.(pjDir)    = saturate(th, opt.ContrastFactorLow(end), opt.ContrastFactorHigh(end), opt.ContrastFactorMethod);
        else
            Pj.LF_th.(pjDir)    = calcThreshold(getFieldVal(Pj.LF,pjDir), opt.ContrastInputCell{1,:});
            Pj.HF_th.(pjDir)    = calcThreshold(getFieldVal(Pj.Alphaval,pjDir)*getFieldVal(Pj.HF,pjDir), opt.ContrastInputCell{1,:});
        end
    end
end

%% Imaging
function [img, Pj] = imaging_internal(LF, HF, AF, opt, Pj, img)
    for pjDir = opt.ProjectionDir
        [x,y]           = axisValues(pjDir, opt);
        Pj              = projection(LF, HF, AF, opt, Pj, pjDir);  
        [f, opt]        = createFig(opt);
        imLF            = applyThresholdAndScaling(Pj.LF.(pjDir), getFieldVal(Pj.LF_th, pjDir), opt);
        switch opt.singleFB
            case true               % Single frequency band  
                imk = imagescCMAP(imLF, opt.CMAP, y, x); 
            case false              % Double frequency band (LF/HF)       
                imHF        = applyThresholdAndScaling(getFieldVal(Pj.Alphaval,pjDir)*Pj.HF.(pjDir), getFieldVal(Pj.HF_th, pjDir), opt);
                if ~isempty(AF)
                    imAF = applyThresholdAndScaling(getFieldVal(Pj.Alphaval,pjDir)*Pj.AF.(pjDir), getFieldVal(Pj.AF_th, pjDir), opt);
                end
                if opt.MergeChannels
                    imMerged = (imLF + imHF)/2;
                    imk = imagescCMAP(imMerged, opt.CMAP, y, x); 
                else
                    imComb= zeros([size(imLF), 3]);
                    rgb_i = 1;
                    for k = opt.ColorChannels
                        switch opt.ColorChannels(rgb_i)
                            case 1
                                imComb(:,:,rgb_i) = imLF;                        
                            case 2
                                imComb(:,:,rgb_i) = imHF;  
                            case 3
                                imComb(:,:,rgb_i) = imAF;  
                        end
                        rgb_i = rgb_i +1;
                    end
                    % imComb      = imfuse(uint8(imLF), uint8(imHF), 'falsecolor','Scaling','none','ColorChannels', opt.ColorChannels); %[1 2 0]
                    imComb = uint8(imComb);
                    imk = image(y,x, imComb);
                end
        end
        setdaspect(opt, pjDir);
        img = setFigProp(opt, img, f, pjDir, imk);
    end
end

%% Figure
function [f , opt] = createFig(opt)
    opt.figureIter = opt.figureIter+1;
    if length(opt.FigureNum) < opt.figureIter
        f = figure();
    elseif ishandle(opt.FigureNum(opt.figureIter))
        set(0, 'CurrentFigure', opt.FigureNum(opt.figureIter));
        f = gcf;
    else
        f = figure(opt.FigureNum(opt.figureIter));
    end
end

function DaspectValue = getDaspectValue(opt, pjDir)
    switch pjDir
    case 'x'
        DaspectValue = opt.DaspectValueX;
    case 'y'
        DaspectValue = opt.DaspectValueY;
    case 'z'
        DaspectValue = opt.DaspectValueZ;
    end

    if isempty(DaspectValue)
        DaspectValue = opt.DaspectValue;
    end

    if ~isempty(opt.AspectXYZ) && (isempty(DaspectValue) || ~isnumeric(DaspectValue))
        axDir = pjDir2axisDir(pjDir, type = "alphabetic");
        DaspectValue = [opt.AspectXYZ([1:axDir-1 axDir+1:end]) 1];
    end

end

function setdaspect(opt, pjDir)
    ax = gca;
    axisDir         = pjDir2axisDir(pjDir);
    DaspectValue   = getDaspectValue(opt, pjDir);
    if ~isnumeric(DaspectValue) && strcmpi(DaspectValue, 'default')
        if axisDir == 1 || strcmp(opt.Unit, 'pixel')
            daspect([1 1 1])
        else
            daspect([2 1 1])
        end
    elseif isstring(DaspectValue) || ischar(DaspectValue)
        axis(ax, DaspectValue)
    elseif iscell(DaspectValue)
        axis(ax, DaspectValue{:})
    elseif isnumeric(DaspectValue) 
        daspect(DaspectValue)
    end
end

function img = setFigProp(opt, img, f, k, imk)
    ax = gca;
    f.Color     = opt.BackgroundColor;
    switch k
    case 'x'
        img.imX = imk;
    case 'y'
        img.imY = imk;
    case 'z'
        img.imZ = imk;
    end
    if ~opt.AxisLabelAndTicks 
        ax.TickDir = "none";
        ax.XTick = [];
        ax.YTick = [];
        f.Visible = opt.Visible;

        return
    end

    ax.FontSize = 14;
    ax.Color    = opt.LabelColor;
    ax.XColor   = opt.LabelColor;
    ax.YColor   = opt.LabelColor;
    ax.ZColor   = opt.LabelColor;
    k           = lower(k);
    fast_slow_str = {'slow scanning','fast scanning'};
    switch k
    case 'x'
        ax.XLabel.String = ['y - ' fast_slow_str{strcmp('y',opt.FastScanningAxis)+1} ' (' opt.Unit ')'];
        ax.YLabel.String = ['z - depth (' opt.Unit ')'];
    case 'y'
        ax.XLabel.String = ['x - ' fast_slow_str{strcmp('x',opt.FastScanningAxis)+1} ' (' opt.Unit ')'];
        ax.YLabel.String = ['z - depth (' opt.Unit ')'];
    case 'z'
        ax.YDir = "normal";
        ax.XLabel.String = ['y - ' fast_slow_str{strcmp('y',opt.FastScanningAxis)+1} ' (' opt.Unit ')'];
        ax.YLabel.String = ['x - ' fast_slow_str{strcmp('x',opt.FastScanningAxis)+1} ' (' opt.Unit ')'];
    end
    set(f, "WindowState", opt.WindowState);
    set(f, "Units", opt.FigureUnits);
    if ~isempty(opt.FigurePosition)
        l = size(opt.FigurePosition,2);
        if size(opt.FigurePosition,1) >= k
            idx = k;
        else
            idx = 1;
        end
        f.Position(1:l) = opt.FigurePosition(idx, 1:l);
    end
    f.Visible = opt.Visible;
end

%% get/set Values
function val = getFieldVal(X, pjDir)
    if isfield(X, pjDir)
        val = X.(pjDir);
    else
        val = X;
    end
end

function [x,y] = axisValues(pjDir, opt)
    switch pjDir
    case 'z'
        % x = [0 opt.ds*size(LF,2)]+opt.ds/2;
        % y = [0 opt.ds*size(LF,3)]+opt.ds/2;
        x = opt.Crop(:,2)'+opt.ds/2;
        y = opt.Crop(:,3)'+opt.ds/2;
    case 'x'
        x = opt.Crop(:,1)'-opt.SurfaceHeight+opt.dz/2;
        y = opt.Crop(:,3)'+opt.ds/2;
    case 'y'
        x = opt.Crop(:,1)'-opt.SurfaceHeight+opt.dz/2;
        y = opt.Crop(:,2)'+opt.ds/2;
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

function mustOnlyContain(a, containChar)
    for ak = a
        if ~any(containChar == ak)
            eidType = 'mustOnlyContain:notOnlyContain';
            msgType = ['Input must only contain following characters: ' containChar];
            throwAsCaller(MException(eidType,msgType))
        end
    end
end
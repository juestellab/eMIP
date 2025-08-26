function [surface, surface_fit, opt] = surfaceDetectionRSOM(R_in, LF, opt)
    % Detects the surface of an RSOM recon by fitting a polynomial of
    % degree 1 to 4. 
    % Input:    - RSOM recon (AF OR HF OR HF/LF)
    %           - Optional: hyper-parameter
    %             Note: input unit is always in pixel
    % Output:   - 2D surface coordinates (height in pixels)
    %           - sfit object (polynomial)
    %           - options used 
    % The algorithm includes:
    %           - Detction of surface points
    %           - Multiple steps to exclude outliers (such as reflection
    %           points, strong absorbing vessles,...)
    %           - the final polynomial fit is performed by the function fit
    %           (2D polynomial of degree 1-4 x 1-3 (fast-scanning x slow-scanning))    
    arguments
        %% User settings
        R_in                                                                                            % HF or full recon        
        LF                                                                  = [];                       % optional LF; if top layer has no low frquency structure it is better to ommit LF; 
        opt.Sensitivity     {mustBeGreaterThanOrEqual(opt.Sensitivity, 1)}  = 1.75;                     % Sensitivity of surface detection; Noisy images -> low sensitivity (e.g., 1.15); 
                                                                                                        % low intensity -> high sensitivity (e.g., 2.5-3.5)     
        opt.SurfaceDepthTop                                                 = 1;                        % Surface height top -> limits the search range from above
        opt.SurfaceDepthBot                                                 = Inf;                      % Surface height bottom -> limits the search range from below
        opt.SelectSurfacePoint            (1,1) logical                     = false;                    % If true, a popup figure appears. There you can select a point closely above the 
                                                                                                        % surface to exclude the negative effect of noise and artifacts above the skin layer.
        %% Devoloper settings
        opt.DispFig char {mustBeMember(opt.DispFig, {'on','off','final'})}  = 'off';                    % if 'on', plots of all iterations of the surface detection will be shwon; if 'final', only the final plot will be shown
        opt.Alpha                          double                           = [0.5  0.9  0.5 0.9];      % alpha qunatiles of the zMIPs for surface point selection (HF low, HF high, LF low, LF high; where low/high means the lower/higher threshold used in the threshold). 
                                                                                                        % Interpretation: choose minimal th_alpha such that 1-alpha <= |pts zMIP >= th_alpha| / |pts zMIP|.
                                                                                                        % 1/Sensitivity * th_alpha * individual thresholdMultiplier = threshold; first point in z direction >= threshold will be selcted as surface point;
        opt.MinMaxAlpha                    double                           = [0.15 0.975 0.5 0.975];   % The final surface selection uses dynamic alphas. MinMaxAlpha gives the range for the HF alpha and the LF alpha, i.e. [minHF_alpha maxHF_alpha minLF_alpha maxLF_alpha]. 
        opt.NoiseRange                                                      = [-20 20];                 % To derive the particular values for the dynamic alphas the background values closely above the surface are used. Noise range defines the lower and uper boundaries of this region (3D volume)
        opt.ThresholdMultiplier            double                           = [.5 1 1 1.25 1 1.25];       % thresholdMultiplier * alpha percentile = threshold; first point in z direction >= threshold will be selcted as surface point; multiplier gives flexibility to the selection process
        opt.RelativeRasterEdgeSize                                          = 0.15;                     % side length (relative) of the raster rectangle (see later explanation)
        opt.AlphaRasterRelaxation         (1,1) double                      = 0.25;                     % alphaRaster = alpha + (1-alpha)*alphaRasterRelaxation (-> in a given raster a local th_alpha is caluclated using the relaxed alpha)
        opt.MinDampingFactor                                                = 0.15;                     % ...the threshold in a raster tile cannot be less than min_damping_factor * (global) surface_threshold (this avoids in regions with no detected points to detect random noise/artifacts)
        opt.ExclusionRatios                     struct                      = struct();                 % see initalization
        opt.SurfaceAbsTolerance           (1,2) double                      = [30 60];                  % (1): choose value below usual height of reflection (above surface) (but not too low); (2) appr. double of first entry
        opt.dzResolutionMultiplyer                                          = 3;                        % resolution in z driection (in um), heights in pixel will be adapted accordingly
        opt.PolyDegreePenaltyFactor                                         = 0.01;                     % penalty factor to penalizes fits using polynomials with high degree
        opt.MaxDegree                     (1,2) double {mustBeInteger(opt.MaxDegree)} = [4 3];          % Maximum degree of surface polynomial (first entry largerer dimension (fast scanning axis))
    end
    %% initalization
    [opt, exc_ratio, exc]               = init(R_in, opt);
    
    %% surface range 00 - manual selection
    [opt, surface_range]                = setSurfaceRange(R_in, LF, opt);

    %% Block 00
    % find surface points 00
    [Global, Raster, exc, RasterLF]     = FBsurfacePtsSelection(0, R_in, LF, surface_range, opt);

    % Set seed
    rng(mod(keyHash(Global),2^32));

    % sfit Iter 00-01 - linear fit 2D polynomial (degree 1) & exclude outliers; 01: linear RANSAC
    [Raster, exc, surface_fit]          = sfitAndExclude(0:1, Global, Raster, exc, exc_ratio, [], opt);

    %% Block 01
    % surface range 01
    surface_range1                      = updateSurfaceRange(1, Global, exc, surface_fit, surface_range, opt);    
    % update surface points 01
    [Global, Raster, exc]               = FBsurfacePtsSelection(1, R_in, LF, surface_range1, opt, RasterLF); 
    % sfit Iter 02-03 - linear fit 2D polynomial (degree 1) & exclude outliers
    [Raster, exc, surface_fit]          = sfitAndExclude(2:3, Global, Raster, exc, exc_ratio, surface_fit, opt);
    % surface range 02
    surface_range                       = updateSurfaceRange(2, Global, exc, surface_fit, surface_range, opt);

    %% Block 02
    % update surface points 02
    [Global, Raster, exc,]              = FBsurfacePtsSelection(2, R_in, LF, surface_range, opt);
    % sfit Iter 04-07 - fit 2D polynomial (degree 1 to 4) & exclude outliers
    [Raster, exc, surface_fit]          = sfitAndExclude(4:7, Global, Raster, exc, exc_ratio, surface_fit, opt);
    % surface range 03
    [surface_range, surface_range_noise]= updateSurfaceRange(3, Global, exc, surface_fit, surface_range, opt);

    %% Block 03
    % adaptive alpha
    alpha_new                           = adaptiveAlpha(R_in, LF, surface_range, surface_range_noise, opt);
    % update surface points 03 - adaptive alpha
    [Global, Raster, exc]               = FBsurfacePtsSelection(3, R_in, LF, surface_range, opt, [], alpha_new);
    % sfit Iter 08:10 - fit 2D polynomial (degree 1 to 4) & exclude outliers; 08: linear RANSAC for correction
    [Raster, exc, surface_fit]          = sfitAndExclude(8:10, Global, Raster, exc, exc_ratio, surface_fit, opt); % fprintf('\n')

    %% correct surface level
    [surface_fit]                       = correctSurfaceLevel(Global, exc, surface_fit, opt);

    %% Surface
    surface = max(1,round(surface_fit(opt.X,opt.Y)));
end
%% Initialization
function [opt, exc_ratio, exc] = init(R_in, opt)
    %% parameters
    opt.dz_resolution_multiplier                = opt.dzResolutionMultiplyer/3;                         % default depth resolution 3um
    surface_range_bdry                          = [1 size(R_in,1)];                                     % can be chosen below the reflection layer to prevent to detect the reflection layer
    opt.surface_depth_top                       = min(max(surface_range_bdry(1), opt.SurfaceDepthTop), surface_range_bdry(2));
    opt.surface_depth_bot                       = min(max(surface_range_bdry(1), opt.SurfaceDepthBot), surface_range_bdry(2));
    AlphaRaster                                 = opt.Alpha + (1-opt.Alpha)*opt.AlphaRasterRelaxation;  % alpha in a given raster (HF low, HF high, LF low, LF high)
    opt.Alpha                                   = [opt.Alpha; AlphaRaster];
    opt.sizeR_in                                = size(R_in);
    exc_fields = {'top_sqrt_secondmoment', 'bot_sqrt_secondmoment', 'top_sigma', 'bot_sigma'};
    for i = 10:-1:1
        switch i
        % SPS (HFhigh, LFhigh, LFlow)    
        % poly11
        case 1
            % remove only points opt.SurfaceAbsTolerance(1) pixels above surface_fit
            exc_ratio(1).top_sqrt_secondmoment  = inf;    
            exc_ratio(1).bot_sqrt_secondmoment  = inf;
            exc_ratio(1).top_sigma              = inf;
            exc_ratio(1).bot_sigma              = inf;
        % RANSAC
        % SPS (HFlow, re-use LFlow) 
        case 2
            % remove points opt.SurfaceAbsTolerance(1) pixels above
            % surface_fit OR 2 sigma below -> remove reflections and noise
            % OR vessels and structures (way) below surface
            exc_ratio(2).top_sqrt_secondmoment  = inf;    
            exc_ratio(2).bot_sqrt_secondmoment  = 0;
            exc_ratio(2).top_sigma              = inf;               
            exc_ratio(2).bot_sigma              = 2;                 
        % poly11
        case 3
            % remove points opt.SurfaceAbsTolerance(2) pixels above
            % surface_fit OR 2 sigma above OR 2 sigma below -> remove reflections and noise
            % OR vessels and structures (way) below surface
            % for linear fits points in each raster are not neccessary for
            % robustness
            exc_ratio(3).top_sqrt_secondmoment  = 0;    
            exc_ratio(3).bot_sqrt_secondmoment  = 0;
            exc_ratio(3).top_sigma              = 2;                
            exc_ratio(3).bot_sigma              = 2;                
        % poly11
        % SPS (HFlow, LFlow, LFlow)
        % poly1-2|1-2
        case 4

            exc_ratio(4).top_sqrt_secondmoment  = 1;    
            exc_ratio(4).bot_sqrt_secondmoment  = 1;
            exc_ratio(4).top_sigma              = 2;
            exc_ratio(4).bot_sigma              = 3;
        % poly1-3|1-3
        case 5
            exc_ratio(5).top_sqrt_secondmoment  = 2;
            exc_ratio(5).bot_sqrt_secondmoment  = 1.5;
            exc_ratio(5).top_sigma              = 2;
            exc_ratio(5).bot_sigma              = 1.5;
        % poly2-3|1-3
        case 6
            exc_ratio(6).top_sqrt_secondmoment  = Inf;
            exc_ratio(6).bot_sqrt_secondmoment  = 1.5;
            exc_ratio(6).top_sigma              = Inf;
            exc_ratio(6).bot_sigma              = 0.75;
        % poly1-3|1-3
        case 7
            exc_ratio(7).top_sqrt_secondmoment  = Inf;
            exc_ratio(7).bot_sqrt_secondmoment  = 2;
            exc_ratio(7).top_sigma              = Inf;
            exc_ratio(7).bot_sigma              = 1;
        % SPS (HFlow, LFlow, LFhigh)
        case 8
            exc_ratio(8).top_sqrt_secondmoment  = 1;
            exc_ratio(8).bot_sqrt_secondmoment  = 1;
            exc_ratio(8).top_sigma              = 3;
            exc_ratio(8).bot_sigma              = 1;
        % RANSAC on resiudal -> update previous poly1-3|1-3 
        case 9
            exc_ratio(9).top_sqrt_secondmoment  = Inf;
            exc_ratio(9).bot_sqrt_secondmoment  = 0.75;
            exc_ratio(9).top_sigma              = Inf;
            exc_ratio(9).bot_sigma              = 0.75;
        % poly1-3|1-3
        case 10
            exc_ratio(10).top_sqrt_secondmoment = Inf;
            exc_ratio(10).bot_sqrt_secondmoment = 0.75;
            exc_ratio(10).top_sigma             = Inf;
            exc_ratio(10).bot_sigma             = 0.75;
        end
        % poly1-3|1-3  

        if size(opt.ExclusionRatios,2) >= i
            for k = exc_fields
                fk = k{1};
                if isfield(opt.ExclusionRatios, fk) && ~isempty(opt.ExclusionRatios(i).(fk))
                    exc_ratio(i).(fk) = opt.ExclusionRatios(i).(fk);
                end
            end
        end
    end
    exc = [];
    % opt.ExclusionRatios = exc_ratio;
    x           = 1:size(R_in,2);
    y           = 1:size(R_in,3);
    [X,Y]       = meshgrid(x,y);
    X           = X';
    Y           = Y';
    opt.X = X;
    opt.Y = Y;
    % Create 4x4 raster of xy with side lengths 15% - 35% - 35% - 15% (with default value raster_sz = 0.15) in x and y
    % ...|.......|.......|...
    % ---x-------x-------x---
    % ...|.......|.......|...
    % ...|.......|.......|...
    % ---x-------x-------x---
    % ...|.......|.......|...
    % ...|.......|.......|...
    % ---x-------x-------x---
    % ...|.......|.......|...

    x_Raster = round(opt.RelativeRasterEdgeSize*size(R_in,2));
    y_Raster = round(opt.RelativeRasterEdgeSize*size(R_in,3));
    x_raster_half = floor((size(R_in,2)+1)/2); 
    y_raster_half = floor((size(R_in,3)+1)/2); 
    raster_vec.x{1} = 1:x_Raster; 
%     raster_vec.x{2} = 1+x_Raster:size(R_in,2)-x_Raster;
    raster_vec.x{2} = 1+x_Raster:x_raster_half;
    raster_vec.x{3} = x_raster_half+1:size(R_in,2)-x_Raster;
    raster_vec.x{4} = size(R_in,2)-x_Raster+1:size(R_in,2);

    raster_vec.y{1} = 1:y_Raster;
%     raster_vec.y{2} = 1+y_Raster:size(R_in,3)-y_Raster; 
    raster_vec.y{2} = 1+y_Raster:y_raster_half;
    raster_vec.y{3} = y_raster_half+1:size(R_in,3)-y_Raster;
    raster_vec.y{4} = size(R_in,3)-y_Raster+1:size(R_in,3);
    opt.raster_vec  = raster_vec;
    opt.SurfaceAbsTolerance = opt.SurfaceAbsTolerance * opt.dz_resolution_multiplier;

    a = size(R_in,2) < size(R_in,3);
    M               = ones(10,2);           N               = ones(10,2); % allows polynomials of degree Mk1 to Mk2 w.r.t. the x-direction in iteration k 
    M(4,      :)    = [1 2];                N(4,      :)    = [1 2];
    M([5,7,9],:)    = repmat([1 3],3,1);    N([5,7,9],:)    = repmat([1 3],3,1);
    M(6,      :)    = [1+(1-a) 3];          N(6,      :)    = [1+a 3];
    M(10,     :)    = [1 opt.MaxDegree(1+a)];          
    N(10,     :)    = [1 opt.MaxDegree(2-a)];

    % M(10,:) = max(M(10,:), opt.MaxDegree(1+a));
    % M(10,:) = max(M(10,:), opt.MaxDegree(1-a));
    opt.polyRangeM = M;                     opt.polyRangeN = N; 
end

%% Surface Points
function [Global, Raster, exc, Raster2] = FBsurfacePtsSelection(iter, R_in, LF, surface_range, opt, RasterOld, alpha) 
    exc     = [];
    Raster2 = [];
    switch iter
        case 0
            thresholdMultiplierHF1  = opt.ThresholdMultiplier(2);           % 2
            thresholdMultiplierLF1  = opt.ThresholdMultiplier(3);           % 3
            thresholdMultiplierLF2  = opt.ThresholdMultiplier(4);           % 4
            alphaHF1                = opt.Alpha(:,2);                       % 2    
            alphaLF1                = opt.Alpha(:,3);                       % 3
            alphaLF2                = opt.Alpha(:,4);                       % 4
        case 1
            thresholdMultiplierHF1  = opt.ThresholdMultiplier(1);           % 1
            alphaHF1                = opt.Alpha(:,1);                       % 1    
        case 2
            thresholdMultiplierHF1  = opt.ThresholdMultiplier(1);           % 1
            thresholdMultiplierLF1  = opt.ThresholdMultiplier(3);           % 3
            alphaHF1                = opt.Alpha(:,1);                       % 1 
            alphaLF1                = opt.Alpha(:,3);                       % 3      
        case 3
            thresholdMultiplierHF1  = opt.ThresholdMultiplier(5);           % 5
            thresholdMultiplierLF1  = opt.ThresholdMultiplier(6);           % 6
            alphaHF1                = alpha(:,1);                           % 1 - dynamic alpha 
            alphaLF1                = alpha(:,end);                         % 2 - dynamic alpha 
    end

    % HF/ ALL frequencies
    [Global, Raster]                = surfacePtsSelection(R_in, surface_range, opt, alphaHF1, thresholdMultiplierHF1);
    
    if isempty(LF)
        return
    end

    % LF
    switch iter
        case 0
            [~, RasterLF1]          = surfacePtsSelection(LF,   surface_range, opt, alphaLF1, thresholdMultiplierLF1);
            [~, RasterLF2]          = surfacePtsSelection(LF,   surface_range, opt, alphaLF2, thresholdMultiplierLF2);
            RasterCell              = {Raster, RasterLF1, RasterLF2};
            Raster2                 = RasterLF1;
        case 1
            RasterCell              = {Raster, RasterOld};
        case 2
            [~, RasterLF1]          = surfacePtsSelection(LF,   surface_range, opt, alphaLF1, thresholdMultiplierLF1);
            RasterCell              = {Raster, RasterLF1, RasterLF1};
        case 3
            [~, RasterLF1]          = surfacePtsSelection(LF,   surface_range, opt, alphaLF1, thresholdMultiplierLF1);
            RasterCell              = {Raster, RasterLF1};
    end
    [Global, Raster]                = appendGlobalAndRaster(RasterCell{:});
end

function [Global, Raster] = surfacePtsSelection(R_in, surface_range, opt, alpha, thMultiplier, surface_threshold)
    if nargin < 5 || isempty(thMultiplier)
        thMultiplier = 1;
    end
    [MIPz, surface_range] = zProjection(R_in, surface_range);    
    if nargin < 6 || isempty(surface_threshold)
        surface_threshold = getThreshold(MIPz, alpha(1,:), 1/thMultiplier*opt.Sensitivity);
    end

    detectedDepth = zeros(size(R_in,2,3));
    Raster = struct();
    RasterSz_x = numel(opt.raster_vec.x);  
    RasterSz_y = numel(opt.raster_vec.y);
    
    for i = RasterSz_x:-1:1
        for j= RasterSz_y:-1:1
            [Raster, detectedDepth] = RasterDetection(R_in, Raster, detectedDepth, MIPz, i, j, surface_range,  alpha, opt, thMultiplier, surface_threshold);
        end
    end
    
    % get indices for all non nan entries (and coordinates)        
    Global.X           = [];
    Global.Y           = [];
    Global.Depth       = [];
    
    for i = 1:RasterSz_x
        for j = 1:RasterSz_y
            detectedDepth_Raster= detectedDepth(opt.raster_vec.x{i},opt.raster_vec.y{j});
            include             = ~isnan(detectedDepth_Raster);
            include             = include(:);
            Raster(i,j).X       = Raster(i,j).X_all(include);
            Raster(i,j).Y       = Raster(i,j).Y_all(include);
            Raster(i,j).Depth   = detectedDepth_Raster(include);
            Global.X            = [Global.X; Raster(i,j).X];
            Global.Y            = [Global.Y; Raster(i,j).Y];
            Global.Depth        = [Global.Depth; Raster(i,j).Depth];
        end
    end
end

function [Raster, detectedDepth] = RasterDetection(R_in, Raster, detectedDepth, MIPz, i, j, surface_range,  alpha, opt, thMultiplier, surface_threshold)
    [X_Raster, Y_Raster]            = meshgrid(opt.raster_vec.x{i}, opt.raster_vec.y{j});
    Raster(i,j).X_all               = X_Raster';
    Raster(i,j).Y_all               = Y_Raster';
    Raster(i,j).exc                 = [];
    MIPz_Raster                     = MIPz(opt.raster_vec.x{i}, opt.raster_vec.y{j});
    th_Raster                       = getThreshold(MIPz_Raster, alpha(2,:), opt.Sensitivity*1/thMultiplier);
    adapted_threshold               = max(th_Raster(end), opt.MinDampingFactor*surface_threshold); % min(max( max(th_Raster(end), opt.MinDampingFactor*surface_threshold), surface_threshold));
    
    % Iterate through every point in Raster
    for ix=opt.raster_vec.x{i}
        for jy=opt.raster_vec.y{j}
            detecedDepth_ij = find(R_in(surface_range(ix,jy,1):surface_range(ix,jy,2),ix,jy,1) > adapted_threshold) + surface_range(ix,jy,1)-1; 
            if isempty(detecedDepth_ij)      
                detecedDepth_ij = nan;
            end
            detectedDepth(ix,jy) = detecedDepth_ij(1);
        end
    end
end

function [MIPz, surface_range] = zProjection(R_in, surface_range)
    if nargin < 2 || isempty(surface_range)
        surface_range           = ones([size(R_in,2,3), 2]);
        surface_range(:,:,2)    = size(R_in,1); 
    end
    surface_range = min(max(surface_range,1), size(R_in,1));
    if max(surface_range(:,:,1),[],'all') == min(surface_range(:,:,1),[],'all') && max(surface_range(:,:,2),[],'all') == min(surface_range(:,:,2),[],'all')
        MIPz                    = squeeze(max(R_in(surface_range(1,1,1):surface_range(1,1,2),:,:), [], 1));
    else
        x = 1:size(R_in,2);
        y = 1:size(R_in,3);
        MIPz = zeros(size(R_in,2,3));
        for ix=x
            for jy=y
                MIPz(ix,jy) = max(R_in(surface_range(ix,jy,1):surface_range(ix,jy,2),ix,jy));
            end
        end
    end
end

function th = getThreshold(MIPz,alpha, Sensitivity)
    MIPz_alpha_percentile   = prctile(MIPz,100*alpha,"all");
    th                      = mean(1/Sensitivity .* MIPz_alpha_percentile);  
end

function [Global, Raster] = appendGlobalAndRaster(Raster, varargin)
    % Combines two or more data sets.
    [N,M] = size(Raster);
    Global.X           = [];
    Global.Y           = [];
    Global.Depth       = [];
    for n = 1:N
        for m = 1:M
            evalX               = @(st) st(n,m).X;
            RX                  = cellfun(evalX, varargin, 'UniformOutput', false);
            evalY               = @(st) st(n,m).Y;
            RY                  = cellfun(evalY, varargin, 'UniformOutput', false);
            evalDepth           = @(st) st(n,m).Depth;
            RDepth              = cellfun(evalDepth, varargin, 'UniformOutput', false);
            Raster(n,m).X       = [Raster(n,m).X; cat(1,RX{:})];
            Raster(n,m).Y       = [Raster(n,m).Y; cat(1,RY{:})];
            Raster(n,m).Depth   = [Raster(n,m).Depth; cat(1,RDepth{:})];
            Global.X            = [Global.X; Raster(n,m).X];
            Global.Y            = [Global.Y; Raster(n,m).Y];
            Global.Depth        = [Global.Depth; Raster(n,m).Depth];
        end
    end
end

%% Adaptive alpha
function alpha_new = adaptiveAlpha(R_in, LF, surface_range, surface_range_noise, opt)
    alpha_new                           = adaptAlpha2Noise(R_in, opt.MinMaxAlpha(1:2), surface_range, surface_range_noise, opt, opt.Sensitivity * opt.ThresholdMultiplier(5));
    % fprintf([num2str(round(alpha_new,3)) ' '])
    if ~isempty(LF)
        alphaLF_new                     = adaptAlpha2Noise(LF, opt.MinMaxAlpha(end-1:end), surface_range, surface_range_noise, opt, opt.Sensitivity * opt.ThresholdMultiplier(6));
        % fprintf([num2str(round(alphaLF_new,3)) ' '])
        alpha_new                       = [alpha_new alphaLF_new];
    end
end

function [alpha_new] = adaptAlpha2Noise(R_in, MinMaxAlpha, surface_range, surface_range_noise, opt, Sensitivity)
    MIPz_noise                          = zProjection(R_in, surface_range_noise);
    MIPz                                = zProjection(R_in, surface_range);
    th_low                              = getThreshold(MIPz, MinMaxAlpha(1), Sensitivity);
    th_high                             = getThreshold(MIPz, MinMaxAlpha(2), Sensitivity);
    th_noise                            = getThreshold(MIPz_noise, 0.95, 1);
    th_new                              = min(max(th_noise, th_low), th_high);
    alpha_new                           = sum(MIPz <= th_new*Sensitivity,'all')/numel(MIPz);
    alphaRaster_new                     = alpha_new + (1-alpha_new)*opt.AlphaRasterRelaxation;
    alpha_new                           = [alpha_new; alphaRaster_new];
end

%% Fit & Exclude
function [Raster, exc, surface_fit] = sfitAndExclude(iterVec, Global, Raster, exc, exc_ratio, surface_fit, opt)
    iter        = iterVec(1);
    % fixed / dynamic surface margin
    switch iter
        case {1,2}
            surface_margin_top              = opt.SurfaceAbsTolerance(1);
        case {3,4,5}
            surface_margin_top              = opt.SurfaceAbsTolerance(2);
        case {6,7,10} % dynamic
            residual_all_points             = surface_fit(Global.X, Global.Y) - Global.Depth;
            residual_interpolation_points   = residual_all_points(~exc);
            surface_margin_top              = prctile(residual_interpolation_points, 80, 'all') + 10*opt.dz_resolution_multiplier;
        case 8
            surface_margin_top              = opt.SurfaceAbsTolerance(2);
            exc                             = zeros(1, length(Global.X));
        case 9
            surface_margin_top              = inf;
    end
    
    % exclude
    if iter > 0
        [exc, Raster]           = exclude_outliers(Raster, exc_ratio(iter), Global, surface_fit, exc, surface_margin_top);
    end

    % plot
    if iter > 0 && iter < 4
        plotSurface(surface_fit, Global, exc, iter-1, opt)
    end

    % fit
    switch iter
        case {0, 2, 3}
            surface_fit         = fit([Global.X,Global.Y], Global.Depth, 'poly11', 'Exclude', exc);
        case {4,5,6,7,9,10}
            surface_fit         = bestPolyfit(Global, exc, opt.polyRangeM(iter,:) , opt.polyRangeN(iter,:), opt, iter);
        case 1
             surface_fit        = applyRansac(Global, exc, 'poly11', 10, opt.SurfaceAbsTolerance(1));
        case 8
            residual_all_points = Global.Depth - surface_fit(Global.X, Global.Y);
            GlobalR.X           = Global.X; 
            GlobalR.Y           = Global.Y; 
            GlobalR.Depth       = residual_all_points;
            surface_fitR        = applyRansac(GlobalR, exc, 'poly11', 10, opt.SurfaceAbsTolerance(1));
            orig_state          = warning;
            warning('off', 'all');
            % fprintf([num2str(round(surface_fitR.p00)) '|' num2str(round(surface_fitR.p01/surface_fit.p01,2)) '|' num2str(round(surface_fitR.p10/surface_fit.p10,2)) ' '])
            surface_fit.p00     = surface_fit.p00 +surface_fitR.p00;
            surface_fit.p01     = surface_fit.p01 +surface_fitR.p01;
            surface_fit.p10     = surface_fit.p10 +surface_fitR.p10;
            warning(orig_state);
    end

    % plot
    switch iter
        case {1}
            plotSurface(surface_fit, Global, exc, 21, opt)
        case {3, 8}
            plotSurface(surface_fit, Global, exc, iter, opt)
    end

    % recursion
    iterVec = iterVec(2:end);
    if ~isempty(iterVec)
        [Raster, exc, surface_fit] = sfitAndExclude(iterVec, Global, Raster, exc, exc_ratio, surface_fit, opt);
    end
end

function [surface_fit] = correctSurfaceLevel(Global, exc, surface_fit, opt)
    residual_all_points                 = surface_fit(Global.X, Global.Y) - Global.Depth;
    residual_interpolation_points       = residual_all_points(~exc);
    offset                              = prctile(residual_interpolation_points,80,'all') + 0.2*std(residual_interpolation_points(residual_interpolation_points>0)) + 3*opt.dz_resolution_multiplier;
    
    orig_state = warning;
    warning('off', 'all');
    surface_fit.p00 = surface_fit.p00 - offset;   
    warning(orig_state);

    if strcmp(opt.DispFig, 'final')
        opt.DispFig = 'on';
    end
    plotSurface(surface_fit, Global, exc, 12, opt)
end

function [surface_fit, usedpoly] = bestPolyfit(Global, exc, M , N, opt, figIter)
    % Runs polynomial fit for each allowed polynomial (e.g. M(1) = 2, N(2) = 3 ->
    % polynomial of the form p(x,y) = p00 + p10*x + p01*y + p20*x^2 +
    % p11*x*y + p02*y^2 + p21*x^2*y + p12*x*y^2 + p03*y^3).
    % Polynomial fit, which minimizes the linear least square to the
    % included surface candidates, will be selcted.
    norm_fit = inf;
    X_c         = Global.X(~exc);
    Y_c         = Global.Y(~exc);
    Depth_c     = Global.Depth(~exc);
    hpartition  = cvpartition(numel(X_c),'Holdout',0.05); 
    idxTrain    = training(hpartition);
    idxTest     = test(hpartition);
    X_Train     = X_c(idxTrain);
    Y_Train     = Y_c(idxTrain);
    Depth_Train = Depth_c(idxTrain);
    X_Test      = X_c(idxTest);
    Y_Test      = Y_c(idxTest);
    Depth_Test  = Depth_c(idxTest);

    orig_state = warning;
    warning('off', 'all');
    for m = M(1):M(2)
        for n = N(1):N(2)
            surface_fit_nm = fit([X_Train,Y_Train], Depth_Train, ['poly' num2str(m) num2str(n)]);
            norm_fit_nm = norm(surface_fit_nm(X_Test, Y_Test) - Depth_Test)*(1+opt.PolyDegreePenaltyFactor)^(norm([m n]-1)^2);
            if norm_fit_nm < norm_fit
                norm_fit = norm_fit_nm;
                usedpoly = ['poly' num2str(m) num2str(n)];
            end
        end
    end
    % fprintf([usedpoly ' '])
    surface_fit = fit([X_c,Y_c], Depth_c, usedpoly);
    warning(orig_state);
    plotSurface(surface_fit, Global, exc, figIter, opt)
end

function surface_fit = applyRansac(Global, exc, usedpoly, sampleSize, maxDistance)
    if isempty(exc)
        Pts                             = [Global.X, Global.Y, Global.Depth];
    else
        Pts                             = [Global.X(~exc),Global.Y(~exc), Global.Depth(~exc)];
    end   
    fitFcn                              = @(x) fit(x(:,1:end-1), x(:,end), usedpoly);
    evalFcn                             = @(model,x) abs((x(:,end)-model(x(:,1:end-1))));
    [surface_fit, inlierIdx, status]    = ransac(Pts,fitFcn,evalFcn, sampleSize, maxDistance);
end

function [exc, Raster] = exclude_outliers(Raster, exc_ratio, Global, surface_fit, exc, varargin)
    % function which detects outlier according to the input parameters
    
    % calculate residual (surface fit-surface candidate depth) of each point globally
    residual        = surface_fit(Global.X, Global.Y) - Global.Depth;
    
    % exclude points above reflection height
    if isempty(exc) && numel(varargin) >= 1 
        exc = residual > varargin{1};     
    end
    if ~isempty(exc)
        residual = residual(~exc);
    end

    % Calculatae global standard deviation of the residual (without
    % excluded points)
    sigma           = std(residual);

    % reset excluded points
    exc = [];
    
    % iterate over each raster (4x4)
    for i = 1:size(Raster,1)
        for j= 1:size(Raster,2)
            % calculate residual (surface fit-surface candidate depth) of
            % each point in raster
            residual        = surface_fit(Raster(i,j).X, Raster(i,j).Y) - Raster(i,j).Depth; 

            % if exclusion vector of raster is not empty, then assign only
            % included (not excluded) candidates to residual_cmp
            if ~isempty(Raster(i,j).exc)
                residual_cmp = residual(~Raster(i,j).exc);
            else
                residual_cmp = residual;
            end

            % calculate sqrt of second moment of residual_cmp in raster, if not
            % availabe -> Inf
            if ~isempty(residual_cmp)
                sqrt_secondmoment = (residual_cmp'*residual_cmp/(numel(residual_cmp)-1))^.5;
            else
                sqrt_secondmoment = inf;
            end
            
            % calculate exclusion threshold (bot/top separately); defined by
            % maximum of ratio_a*sqrt_secondmoment and ratio_b*sigma
            % (Why? Globally, sigma can be low, but in a certain raster the
            % second_moment can be high -> too many candidates 
            % could be excluded by using only sigma)
            exc_bot = residual < -max(exc_ratio.bot_sqrt_secondmoment * sqrt_secondmoment, exc_ratio.bot_sigma * sigma);
            exc_top = residual >  max(exc_ratio.top_sqrt_secondmoment * sqrt_secondmoment, exc_ratio.top_sigma * sigma);
            
            % if reflection height is given, exclude surface height with
            % resiudal bigger than opt.SurfaceAbsTolerance
            if numel(varargin) >= 1
                exc_reflection = residual > varargin{1};
                exc_top = exc_top | exc_reflection;
            end
            
            % assign exclusion
            Raster(i,j).exc = exc_bot | exc_top;
            exc             = [exc; Raster(i,j).exc];
        end
    end
end

%% SurfaceRange
function [opt, surface_range] = setSurfaceRange(R_in, LF, opt)
    if opt.SelectSurfacePoint
        fastDir = {'y','x'}; fastDir = fastDir{(size(R_in,2)>size(R_in,3))+1};
        if isempty(LF)
            imshowRSOM(R_in, "CMAP", 2, "unit", "pixel", "daspectValue", [1 4 1], "surfaceHeightInPixel", 0, "ProjectionDir", fastDir);
        else
            imshowRSOM(LF, R_in, "unit", "pixel", "daspectValue", [1 4 1], "surfaceHeightInPixel", 0, "ProjectionDir", fastDir);
        end
        s = gcf;
        s.Name = "Select point closely above surface!";
        h = drawcrosshair;
        opt.surface_depth_top = min(max(1,round(h.Position(2))),size(R_in,1)-3);
        close(s)      
    end

    surface_range                               = zeros([size(R_in,2,3),2]);
    surface_range(:,:,1)                        = opt.surface_depth_top;
    surface_range(:,:,2)                        = opt.surface_depth_bot;
end

function [surface_range, surface_range_noise] = updateSurfaceRange(iter, Global, exc, surface_fit, surface_range, opt)
    surface_range_noise                         = [];
    surface                                     = max(1,surface_fit(opt.X,opt.Y));
    switch iter
        case 1
            surface_range(:,:,1)               = round(max(surface_range(:,:,1), surface-opt.SurfaceAbsTolerance(2) * opt.dz_resolution_multiplier));
        case 2
            residual_all_points                 = surface_fit(Global.X, Global.Y) - Global.Depth;
            residual_interpolation_points       = residual_all_points(~exc);
            offset                              = prctile(residual_interpolation_points, 90, 'all');
            surface_range(:,:,1)                = round(max(opt.surface_depth_top, surface-opt.SurfaceAbsTolerance(1)-offset));
        case 3
            surface                             = max(1,surface_fit(opt.X,opt.Y));
            surface_range(:,:,1)                = round(min(opt.surface_depth_bot - 1,max(opt.surface_depth_top, surface - opt.SurfaceAbsTolerance(2) * opt.dz_resolution_multiplier)));
            surface_range(:,:,2)                = round(min(opt.surface_depth_bot, surface + opt.SurfaceAbsTolerance(1) * opt.dz_resolution_multiplier));
            surface_range_noise                 = surface_range + opt.NoiseRange(1); 
            surface_range_noise(:,:,2)          = surface_range(:,:,1)+opt.NoiseRange(2);
            surface_range_noise                 = max(1, min(opt.sizeR_in(1), surface_range_noise));
    end
end

%% Visualization
function plotSurface(surface_fit, Global, exc, figIter, opt)
    if strcmp(opt.DispFig, 'on')
        figure(200+figIter)
        plot(surface_fit, [Global.X, Global.Y], Global.Depth, 'Exclude', exc)
        set(gca, ZDir = 'reverse')
        if figIter == 12
            figure(213)
            plot(surface_fit, [Global.X(~exc), Global.Y(~exc)], Global.Depth(~exc))
            set(gca, ZDir = 'reverse')
        end
    end
end

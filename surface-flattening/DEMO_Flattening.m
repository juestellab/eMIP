%% DEMO flattening
% [filename,folder_path ] = uigetfile('Z:\RSOM_Data\RSOM_backup_Suhanyaa\DiabetesRecon\HF_LF','Select RSOM reconstruction.');
[filename,folder_path ] = uigetfile('C:','Select RSOM reconstruction.');

filenameLF = replace(filename, 'HF', 'LF');
filenameHF = replace(filename, 'LF', 'HF');
HF = load([folder_path, filenameHF]);
LF = load([folder_path, filenameLF]);

%% Standard Flattening
 [surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R, "dispFig", "on");
 imshowRSOM(LF,HF, "Surface", surface_fit)

%% IF only low frequency structure is visible (e.g., psoriasis data sets) 
 [surface, surface_fit] = surfaceDetectionRSOM(HF.R,  "dispFig", "final");
 imshowRSOM(LF,HF, "Surface", surface_fit);


 %% Change Sensitivity if surface is not detected correctly 
 % (Noisy images -> low sensitivity (e.g., 1.15));
 [surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R, "dispFig", "final", "Sensitivity", 1.15);
 imshowRSOM(LF,HF,  "Surface", surface_fit);

 %% Select surface point
 % the surface height can be selected manually to improve results

 [surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R, "dispFig", "final", "SelectSurfacePoint", true);
 imshowRSOM(LF,HF,  "Surface", surface_fit);


%% Create new corrected 3d-Volume
HF_R_flat = surfaceCorrectionRSOM(HF.R, surface);

%% change surfcae position (in pixel)
HF_R_flat = surfaceCorrectionRSOM(HF.R, surface, "SurfacePos", "top");
imshowRSOM(HF,  "SurfaceHeightInmm", 0, "CMAP","gray", "ProjectionDir","x");
imshowRSOM(HF,  "Surface", surface_fit, "SurfaceHeightInmm", 0, "SurfaceVarargin", {"SurfacePos", "top"}, "CMAP","gray", "ProjectionDir","x");


%% correct and make sure all data is preserved (with zero paadding)
HF_R_flat = surfaceCorrectionRSOM(HF.R, surface, "ImageSize", "max", "SurfacePos", "bottom");
imshowRSOM(HF,  "Surface", surface_fit,"SurfaceHeightInmm", 0,  "SurfaceVarargin", {"ImageSize", "max", "SurfacePos", "bottom"}, "CMAP","gray", "ProjectionDir","x");



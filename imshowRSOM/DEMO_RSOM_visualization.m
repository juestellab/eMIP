%% DEMO for visualizing RSOM data
% Author: Manuel Gehmeyr
% Email: manuel.gehmeyr@helmholtz-munich.de
%% Load RSOM reconstruction (LF/HF)
[filename,folder_path ] = uigetfile('','Select RSOM reconstruction.');
filenameLF = replace(filename, 'HF', 'LF');
filenameHF = replace(filename, 'LF', 'HF');
HF = load([folder_path, filenameHF]);
LF = load([folder_path, filenameLF]);

%% default visualzation
% Maximum intensity projections (MIPs) along the x-, y- and z-axis for low 
% frequency (LF) and high frequency (HF) reconstruction. The contrast is 
% calculated for LF and HF separately by clipping negative values to 0 and 
% capping valuesat 1.25 times the 95th percentile value of the 
% z-MIP. The corresponding MIPs are then fused into three rgb-images. LF is
% assigned linearly scaled to the red channel and HF to the green channel. 
% Output: 
% -img: structure with image objects, e.g.; img.imZ.CData shows the
% 3D-RGB data of the shown z-MIP.
% - Pj: projection data and infromation
% - opt: options to create the plots

[img, Pj, opt] = imshowRSOM(LF, HF);

%% Adapt Contrast
img = imshowRSOM(LF, HF, 'ContrastAlphaHigh', 0.9, 'ContrastAlphaLow',0.1, 'ContrastFactorHigh', 1, 'ContrastFactorLow', -0.05, 'ContrastFactorMethod','absolute');

%% Assign Figure Number & Define projection Directions
img = imshowRSOM(LF, HF, 'FigureNum', [7, 8], 'ProjectionDir', 'zx');

%% Merge into monochromatic channel
img = imshowRSOM(LF, HF, 'MergeChannels', true, 'CMAP', 'gray',  'FigureNum', [31, 32, 33]);

%% Single channel visualization
img = imshowRSOM(LF, [] , 'CMAP', 'r',  'FigureNum', [31, 32, 33]);

%% Plain image, white background
img = imshowRSOM(LF, HF, 'ProjectionDir', 'x', 'FigureNum',1, 'AxisLabelAndTicks', false, 'BackgroundColor','w');

%% Units and aspect ratio
img = imshowRSOM(LF, HF, 'ProjectionDir', 'zx', 'DaspectValueX', [1 1 1], 'DaspectValueZ', [2 1 1] ,"Unit", "pixel");   % see matlab documentation daspect
img = imshowRSOM(LF, HF, 'ProjectionDir', 'zx', 'DaspectValueX', [1 1 1], 'DaspectValueZ', [2 1 1] ,"Unit", "um");      % see matlab documentation daspect

%% scaling
img = imshowRSOM(HF, [], 'ContrastScaling', 'sqrt',                                                         'ProjectionDir', 'x', 'CMAP', 'gray');
img = imshowRSOM(HF, [], 'ContrastScaling', 'log',   'ContrastFactorHigh', 10,   'ContrastFactorLow', 0.01, 'ProjectionDir', 'x', 'CMAP', 'gray');
img = imshowRSOM(HF, [], 'ContrastScaling', 'power', 'ContrastScalingPower', 2,  'ContrastFactorLow', -.5,  'ProjectionDir', 'x', 'CMAP', 'gray');
img = imshowRSOM(HF, [], 'ContrastScaling', 'linear',                                                       'ProjectionDir', 'x', 'CMAP', 'gray');

%% croping
img = imshowRSOM(LF,HF, Crop= [450; 650], Unit="um")

%% Fequency band equalization
% [1] J. Aguirre, M. Schwarz, N. Garzorz, M. Omar, A. Buehler, K. Eyerich, 
% and V. Ntziachristos, “Precision assessment of label-free psoriasis 
% biomarkers with ultra-broadband optoacoustic mesoscopy,” Nat. Biomed. 
% Eng., vol. 1, no. 5, May 2017.

img = imshowRSOM(LF, HF, FBE = true, figureNum = [34,35,36]);
% img = imshowRSOM(LF, HF, contrast_low=0.025, figureNum = [34, 35, 36])
imshowRSOM(LF, HF, figureNum = [34,35,36], FBE=true, ContrastAlphaLow= 0, ContrastAlphaHigh= 1, ContrastFactorHigh= 0.35, ContrastFactorLow= 0.06, FBEoptions = 'Pji', BackgroundColor='w', AxisLabelAndTicks=false);

%% Visualize flattened data
[surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R);
[img, Pj, opt] = imshowRSOM(LF, HF, 'Surface', surface_fit, 'SurfaceHeightInmm', 0.5, 'FigureNum', [1 2 3]);


%% Save figures
example_path = './RSOM_example_MIP';
exportgraphicsRSOM(example_path, [1 2 3], {'.png', '.pdf'})

%% Save figures without showing them - usful for batch processes to run in the background
example_path2 = './RSOM_example2_MIP';
imshowRSOM(LF, HF, 'FigureNum', [4 5 6], 'Visible', 'off');
exportgraphicsRSOM(example_path2, [4 5 6])
close all


# eMIP — Enhanced Maximum Intensity Projection for Optoacoustic Imaging (MATLAB)

**eMIP** is a MATLAB toolbox for improving the fidelity and interpretability of optoacoustic images, with a focus on **Raster-Scan Optoacoustic Mesoscopy (RSOM)**.  

Conventional maximum intensity projections (MIPs) often suffer from distortion due to skin curvature, reflections, and depth compression. eMIP introduces methods to **detect and flatten the skin surface**, remove outliers, and generate **enhanced depth-resolved projections** that better represent vascular and structural features.

---

## Features
- Automatic **skin surface detection** with polynomial fitting and RANSAC outlier removal  
- **Surface flattening** and corrected volume generation  
- Flexible **RSOM visualization**: MIPs, multi-channel fusion, contrast scaling, cropping, frequency band equalization (FBE), and inline flattening  
- **Demo scripts** for quick evaluation and reproducibility  

---

## Requirements
- **MATLAB R2023a or higher**  
  (Surface detection uses RANSAC, which changed implementation between releases.)  
- **Toolboxes**:
  - Curve Fitting Toolbox  
  - Machine Learning Toolbox **or** Statistics and Machine Learning Toolbox  

---

## Installation

Clone the repository and add it to your MATLAB path:

```bash
git clone https://github.com/juestellab/eMIP.git
```

In MATLAB:

```matlab
addpath(genpath('eMIP'))
```

---

## Expected Input Data
The demos assume **paired HF/LF reconstructions** in `.mat` files with fields:
- `HF.R` — 3D RSOM volume (z × x × y) for **high frequency**
- `LF.R` — 3D RSOM volume (z × x × y) for **low frequency**
  
---

## Demo Scripts

- **`DEMO_flattening.m`**  
  Demonstrates automatic surface detection and flattening of RSOM data. Includes examples for noisy images, psoriasis datasets (HF-only), sensitivity adjustment, manual point selection, and volume correction with/without padding.

- **`DEMO_RSOM_visualization.m`**  
  Comprehensive visualization workflow for RSOM reconstructions. Includes:
  - MIPs along x/y/z  
  - Contrast scaling (linear, log, sqrt, power)  
  - LF/HF fusion (RGB) or monochrome colormaps  
  - Aspect ratio, cropping, background options  
  - Frequency Band Equalization (FBE)  
  - Inline flattening and figure exports  

Run them directly in MATLAB to see the pipeline in action.

---

## Usage Example

```matlab
% Input: RSOM reconstruction (HF and LF)
[filename, folder] = uigetfile('', 'Select RSOM reconstruction');
HF = load(fullfile(folder, strrep(filename,'LF','HF')));
LF = load(fullfile(folder, strrep(filename,'HF','LF')));

% 1) Detect surface
[surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R, "DispFig","final");

% 2) Correct volume
HF_flat = surfaceCorrectionRSOM(HF.R, surface, "ImageSize","max","SurfacePos","top");

% 3) Visualize flattened data
imshowRSOM(LF, HF, "Surface", surface_fit, "Unit","um", "CMAP","gray");
```

---

## Function Reference

The eMIP toolbox is built around three core functions:

### `surfaceDetectionRSOM`

```matlab
[surface, surface_fit, opt] = surfaceDetectionRSOM(R_in, LF, opt)
```

Detects the skin surface in an RSOM reconstruction by fitting a 2D polynomial surface.  
Includes multiple outlier exclusion steps (reflections, vessels, artefacts).

- **Inputs:**  
  - `R_in` — 3D RSOM volume (HF or full)  
  - `LF` — *(optional)* LF volume (omit if top layer invisible)  
  - `opt` — name–value args  

- **Outputs:**  
  - `surface` — 2D surface coordinates (pixels)  
  - `surface_fit` — polynomial fit (`sfit`)  
  - `opt` — resolved options  

- **Key options:**  
  - `Sensitivity` (**1.75**) — lower for noisy, higher for low intensity  
  - `SelectSurfacePoint` (**false**) — interactively pick a point above surface  

- **Examples:**
```matlab
[surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R, "DispFig","final");
[surface, surface_fit] = surfaceDetectionRSOM(HF.R, "DispFig","final");
[surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R, "Sensitivity",1.15, "SelectSurfacePoint",true);
```

---

### `surfaceCorrectionRSOM`

```matlab
R_flat = surfaceCorrectionRSOM(R_in, surface, opts)
```

Flattens an RSOM volume to the detected surface, optionally with padding to preserve all voxels.

- **Inputs:**  
  - `R_in` — 3D RSOM volume  
  - `surface` — 2D surface coordinates  
  - `opts` — name–value args  

- **Output:**  
  - `R_flat` — flattened volume  

- **Key options:**  
  - `SurfacePos` (**'bottom'**) | `'top'`  
  - `ImageSize` (**'input'**) | `'max'` (preserve all voxels via padding)  

- **Examples:**
```matlab
HF_flat = surfaceCorrectionRSOM(HF.R, surface);
HF_flat = surfaceCorrectionRSOM(HF.R, surface, "ImageSize","max","SurfacePos","top");
imshowRSOM(HF, "Surface", surface_fit, "SurfaceVarargin",{"ImageSize","max","SurfacePos","bottom"});
```

---

### `imshowRSOM`

```matlab
[img, Pj, opt, R_C] = imshowRSOM(LF, HF, opts)
```

Visualize RSOM reconstructions with MIPs, color fusion, scaling, cropping, flattening, and exports.  

- **Inputs:**  
  - `LF`, `HF` — structs with `.R` volumes (can be `[]`)  
  - `opts` — name–value args  

- **Outputs:**  
  - `img` — figure/image handles (`img.imZ.CData`)  
  - `Pj` — projection data  
  - `opt` — resolved options  
  - `R_C` — processed volume  

- **Common options:**  
  - **Layout**: `FigureNum`, `ProjectionDir`, `Visible`, `BackgroundColor`, `AxisLabelAndTicks`  
  - **Units**: `Unit`, `DaspectValueX`, `DaspectValueZ`  
  - **Contrast**: `ContrastScaling`, `ContrastScalingPower`, `ContrastAlphaHigh/Low`, `ContrastFactorHigh/Low`  
  - **Channels**: `MergeChannels`, `CMAP`  
  - **Cropping**: `Crop=[zmin;zmax]`  
  - **FBE**: `FBE`, `FBEoptions`  
  - **Flattening**: `Surface`, `SurfaceHeightInmm`, `SurfaceVarargin`  

- **Examples:**
```matlab
imshowRSOM(LF, HF);
imshowRSOM(HF, [], 'ContrastScaling','log','ContrastFactorHigh',10,'ContrastFactorLow',0.01);
imshowRSOM(LF, HF, "Surface",surface_fit, "SurfaceHeightInmm",0.5, "Crop",[450;650], "Unit","um","CMAP","gray");
```

---

⚡ **Typical Workflow**
```matlab
[surface, surface_fit] = surfaceDetectionRSOM(HF.R, LF.R);
HF_flat = surfaceCorrectionRSOM(HF.R, surface, "ImageSize","max");
imshowRSOM(LF, HF, "Surface",surface_fit,"Unit","um","CMAP","gray");
```

---

## Repository Structure
```
eMIP/
├─ surface-flattening/
│  ├─ DEMO_flattening.m
│  └─ (flattening functions)
├─ visualization/
│  ├─ DEMO_RSOM_visualization.m
│  └─ (imshowRSOM, exportgraphicsRSOM, helpers)
├─ docs/
└─ LICENSE
```

---

## Authors & Manuscript
This toolbox accompanies the manuscript:  

**"Enhanced Maximum Intensity Projection (eMIP) for Improving the Fidelity of Optoacoustic Images"**  
*Manuel Gehmeyr\*, María Begona Rojas López\*, Suhanyaa Nitkunanantharajah, Hubert Preißl, Andreas Vosseler, Reiner Jumpertz von Schwartzenberg, Andreas L. Birkenfeld, Nikoletta Katsouli, Nikolina-Alexia Fasoula, Angelos Karlas, Michael Kallmayer, Anette-Gabriele Ziegler, Dominik Jüstel, Vasilis Ntziachristos‡*  

\* Equal contribution  
‡ Corresponding author: bioimaging.translatum@tum.de  

📄 Published in npj Imaging (2025) 3, 49
https://doi.org/10.1038/s44303-025-00112-z

---

## Contributing
- Open an [issue](https://github.com/juestellab/eMIP/issues) to report bugs or suggest features.  
- Submit a pull request with improvements or additional modules.  

---

## License
This project is released under the [MIT License](LICENSE).  
You are free to use, modify, and distribute this software with attribution.

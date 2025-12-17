# 4Pi-BRAINSPOT Toolbox

**Interferometric Ultra-High Resolution 3D Imaging through Brain Sections**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b-blue.svg)](https://www.mathworks.com/)
[![CUDA](https://img.shields.io/badge/CUDA-7.5+-green.svg)](https://developer.nvidia.com/cuda-downloads)

## 📖 General Information

[cite_start]The **4Pi-BRAINSPOT toolbox** is a software package developed for **in situ point spread function (PSF) retrieval** with 4Pi single-molecule switching nanoscopy (4Pi-SMSN)[cite: 13, 14].

[cite_start]It captures *in situ* interferometric information directly from acquired single-molecule datasets, allowing for the generation of an *in situ* 4Pi-PSF model that minimizes data-model disparity[cite: 15]. [cite_start]The toolbox features a user-friendly interface that covers the entire 4Pi reconstruction workflow[cite: 16].

## ✨ Key Features

* [cite_start]**Complete Workflow:** Handles channel alignment, PSF segmentation, *in situ* retrieval, dynamic model updates, and super-resolution reconstruction[cite: 16].
* [cite_start]**In Situ Modeling:** Constructs 3D PSF models directly from experimental data[cite: 39].
* [cite_start]**Dynamic Correction:** Estimates cavity-phase-induced and objective-misaligned interferometric aberrations[cite: 40].
* [cite_start]**GPU Acceleration:** Supports pupil-based 3D localization using Cubic Spline GPU implementation[cite: 16, 509].
* [cite_start]**Post-Processing:** Includes 3D drift correction, volume alignment, and data visualization[cite: 41, 42].
* [cite_start]**Camera Support:** Supports both EMCCD and sCMOS camera modes (with calibration)[cite: 143].

## 💻 System Requirements

To ensure the software runs correctly, please verify the following environment:

* [cite_start]**Operating System:** Windows 7 or later (64-bit)[cite: 19].
* [cite_start]**Software:** MATLAB R2019b (64-bit)[cite: 20].
* [cite_start]**Hardware:** CUDA 7.5 compatible graphics driver (required for GPU-based localization)[cite: 21].

## 🚀 Installation

1.  **Unzip the Package:**
    Unzip the `4Pi-BRAINSPOT toolbox.zip` file. [cite_start]You will see three folders: `4Pi-BRAINSPOT toolbox`, `Support`, and `Data_4Pi`[cite: 23].

2.  **Configure Paths:**
    [cite_start]Open the `main.m` file located in the `4Pi-BRAINSPOT toolbox` folder[cite: 24]. [cite_start]Check the `support path` variable to ensure it correctly points to the `Support` folder[cite: 25].

3.  **Launch:**
    [cite_start]Run the `main.m` file in MATLAB to start the GUI[cite: 26].

## 🛠️ Step-by-Step Guide

[cite_start]The interface is divided into eight modules[cite: 33]:

### [cite_start]1. Setup [cite: 121]
Configure the workspace and general parameters, including:
* [cite_start]Pixel size and Numerical Aperture (NA)[cite: 136, 140].
* [cite_start]Refractive indices (immersion and sample media)[cite: 137, 138].
* [cite_start]Camera settings (Gain, Offset, and optional sCMOS calibration)[cite: 141, 143].

### [cite_start]2. Data Import [cite: 177]
[cite_start]Import single-molecule interferometric datasets from the 4 detection channels (p1, s2, p2, s1)[cite: 178].
* [cite_start]*Option:* Enable background subtraction using temporal median filtering[cite: 183].

### [cite_start]3. Channel Alignment [cite: 211]
Align image stacks from the 4 detection channels to the same region of interest.
* [cite_start]Calculate affine transformations or import an existing calibration file (`tform`)[cite: 215, 216].

### [cite_start]4. 4Pi-PSF Segmentation [cite: 241]
Crop pairs of sub-regions containing single molecules to construct a PSF library.
* [cite_start]Adjust `Box size`, `Dist thresh`, and `Seg thresh` to ensure isolated molecules are selected[cite: 243, 245, 248].

### [cite_start]5. 4Pi In Situ PSF Retrieval [cite: 291]
Generate an *in situ* 4Pi-PSF model directly from the segmented sub-regions.
* [cite_start]This module uses iterative coherent 4Pi phase retrieval to estimate pupil functions[cite: 293].
* [cite_start]Outputs: Retrieved 4Pi-PSFs, pupil magnitude/phase, and Zernike coefficients[cite: 323, 324, 325].

### [cite_start]6. Dynamic Model Update [cite: 408]
Estimate time-varying interferometric aberrations caused by cavity phase shifts and objective misalignment.

### [cite_start]7. 4Pi Localization [cite: 474]
Reconstruct the 3D super-resolution image.
* [cite_start]**Note:** This module requires a CUDA-compatible GPU[cite: 524].
* [cite_start]Includes segmentation, rejection (based on photon count, LLR, uncertainty), 3D drift correction, and volume alignment[cite: 475, 505].

### [cite_start]8. Display [cite: 604]
Generate x-y views of the reconstructed volume where molecules are color-coded by their axial (z) position.

## 📂 Demo Dataset

[cite_start]A demonstration dataset is provided in the `Data_4Pi` folder[cite: 639]:

* [cite_start]`rawData_4Pi.mat`: Sample 4Pi single-molecule dataset[cite: 640].
* [cite_start]`config_4Pi.mat`: General setting parameters[cite: 641].
* [cite_start]`tform_all.mat`: Alignment calibration file[cite: 642].
* [cite_start]`SCMOS_calibration_4Pi.mat`: sCMOS calibration parameters[cite: 643].
* [cite_start]`recon4Pi.mat`: Sample reconstruction results[cite: 644].

## 👥 Authors

[cite_start]This toolbox is accompanying software for the manuscript **"Interferometric Ultra-High Resolution 3D Imaging through Brain Sections"**[cite: 2].

* [cite_start]**Affiliations:** Purdue University & Beijing Institute of Technology[cite: 4].
* [cite_start]**Correspondence:** Fan Xu, Alexander A. Chubykin, Fang Huang[cite: 9].

[cite_start]For further updates, please check the associated Github repository[cite: 29].
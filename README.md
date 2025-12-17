# 4Pi-BRAINSPOT Toolbox

**Interferometric Ultra-High Resolution 3D Imaging through Brain Sections**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b-blue.svg)](https://www.mathworks.com/)
[![CUDA](https://img.shields.io/badge/CUDA-7.5+-green.svg)](https://developer.nvidia.com/cuda-downloads)

## 📖 General Information

The **4Pi-BRAINSPOT toolbox** is a software package developed for **in situ point spread function (PSF) retrieval** with 4Pi single-molecule switching nanoscopy (4Pi-SMSN).

It captures *in situ* interferometric information directly from acquired single-molecule datasets, allowing for the generation of an *in situ* 4Pi-PSF model that minimizes data-model disparity. The toolbox features a user-friendly interface that covers the entire 4Pi reconstruction workflow.

## ✨ Key Features

* **Complete Workflow:** 4Pi channel alignment, PSF segmentation, *in situ* retrieval, dynamic model updates, and super-resolution reconstruction.
* **In Situ Modeling:** Constructs 3D PSF models directly from experimental data.
* **Dynamic Correction:** Estimates cavity-phase-induced and objective-misaligned interferometric aberrations.
* **GPU Acceleration:** Supports pupil-based 3D localization using Cubic Spline GPU implementation.
* **Post-Processing:** Includes 3D drift correction, volume alignment, and data visualization.
* **Camera Support:** Supports both EMCCD and sCMOS camera modes (with calibration).

## 💻 System Requirements

To ensure the software runs correctly, please verify the following environment:

* **Operating System:** Windows 7 or later (64-bit).
* **Software:** MATLAB R2019b (64-bit).
* **Hardware:** CUDA 7.5 compatible graphics driver (required for GPU-based localization).

## 🚀 Installation

1.  **Unzip the Package:**
    Unzip the `4Pi-BRAINSPOT toolbox.zip` file. You will see three folders: `4Pi-BRAINSPOT toolbox`, `Support`, and `Data_4Pi`.

2.  **Configure Paths:**
    Open the `main.m` file located in the `4Pi-BRAINSPOT toolbox` folder. Check the `support path` variable to ensure it correctly points to the `Support` folder.

3.  **Launch:**
    Run the `main.m` file in MATLAB to start the GUI.

## 🛠️ Step-by-Step Guide

The interface is divided into eight modules:

### 1. Setup
Configure the workspace and general parameters, including:
* Pixel size and Numerical Aperture (NA).
* Refractive indices (immersion and sample media).
* Camera settings (Gain, Offset, and optional sCMOS calibration).

### 2. Data Import
Import single-molecule interferometric datasets from the 4 detection channels (p1, s2, p2, s1).
* *Option:* Enable background subtraction using temporal median filtering.

### 3. Channel Alignment
Align image stacks from the 4 detection channels to the same region of interest.
* Calculate affine transformations or import an existing calibration file (`tform`).

### 4. 4Pi-PSF Segmentation
Crop pairs of sub-regions containing single molecules to construct a PSF library.
* Adjust `Box size`, `Dist thresh`, and `Seg thresh` to ensure isolated molecules are selected.

### 5. 4Pi In Situ PSF Retrieval
Generate an *in situ* 4Pi-PSF model directly from the segmented sub-regions.
* This module uses iterative coherent 4Pi phase retrieval to estimate pupil functions.
* Outputs: Retrieved 4Pi-PSFs, pupil magnitude/phase, and Zernike coefficients.

### 6. Dynamic Model Update
Estimate time-varying interferometric aberrations caused by cavity phase shifts and objective misalignment.

### 7. 4Pi Localization
Reconstruct the 3D super-resolution image.
* **Note:** This module requires a CUDA-compatible GPU.
* Includes segmentation, rejection (based on photon count, LLR, uncertainty), 3D drift correction, and volume alignment.

### 8. Display
Generate x-y views of the reconstructed volume where molecules are color-coded by their axial (z) position.

## 📂 Demo Dataset

A demonstration dataset is provided in the `Data_4Pi` folder:
* `rawData_4Pi.mat`: Sample 4Pi single-molecule dataset.
* `config_4Pi.mat`: General setting parameters.
* `tform_all.mat`: Alignment calibration file.
* `SCMOS_calibration_4Pi.mat`: sCMOS calibration parameters.
* `recon4Pi.mat`: Sample reconstruction results.

## 👥 Authors

This toolbox is accompanying software for the manuscript **"Interferometric Ultra-High Resolution 3D Imaging through Brain Sections"**.

* **Affiliations:** Purdue University & Beijing Institute of Technology.
* **Correspondence:** Fan Xu, Alexander A. Chubykin, Fang Huang.

For further updates, please check the associated Github repository.